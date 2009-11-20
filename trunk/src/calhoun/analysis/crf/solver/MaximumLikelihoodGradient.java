/**
 * 
 */
package calhoun.analysis.crf.solver;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.CRFObjectiveFunctionGradient;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.analysis.crf.solver.CacheProcessor.FeatureEvaluation;
import calhoun.analysis.crf.solver.CacheProcessor.SolverSetup;
import calhoun.analysis.crf.solver.check.AllSparseLengthCacheProcessor;
import calhoun.util.Assert;
import calhoun.util.ColtUtil;

/** computes the likelihood of the true path for a Markov CRF.  The likelihood is normalized to a per label likelihood so
 * that likelihood of different length paths can be meaningfully compared and a single set of optimization tolerances can be used.
 * Must be configured with a {@link CacheProcessor}. */
public class MaximumLikelihoodGradient implements CRFObjectiveFunctionGradient {
	private static final Log log = LogFactory.getLog(MaximumLikelihoodGradient.class);
	boolean debug = log.isDebugEnabled();

	CacheProcessor cacheProcessor = new AllSparseLengthCacheProcessor();
	SolverSetup modelInfo;
	FeatureEvaluation[] evals;
	boolean[] invalidTransitions;
	
	/// Cached value of the Mi matrix for all of the features present at every position
	/// Mi is stored as a sparse matrix
	int miLength;
	double[] mi;

	int iter = 0;

	double[] prevAlpha;
	double[] alpha;
	double[][] betas;
	double[] betaNorms;
	double[] expects;

	// We publish feature sums 
	private double[] featureSums;
	
	/** gets the cache processor used to access feature evaluations
	 * @return the configured cache processor
	 */
	public CacheProcessor getCacheProcessor() {
		return cacheProcessor;
	}

	/** sets the cache processor used to access feature evaluations
	 * @param cacheProcessor the cache processor to use
	 */
	public void setCacheProcessor(CacheProcessor cacheProcessor) {
		this.cacheProcessor = cacheProcessor;
	}

	public void setTrainingData(ModelManager fm, List<? extends TrainingSequence<?>> data) {
		cacheProcessor.setTrainingData(fm, data);
		modelInfo = cacheProcessor.getSolverSetup();
		evals = cacheProcessor.getFeatureEvaluations();
		invalidTransitions = cacheProcessor.getInvalidTransitions();
		
		miLength = modelInfo.nTransitions;
		mi = new double[miLength];

		expects = new double[modelInfo.nFeatures];
		prevAlpha = new double[modelInfo.nStates];
		alpha = new double[modelInfo.nStates];

		// Initialize betas (for use later, in the gradient computation)
		betas = new double[modelInfo.longestSeq][modelInfo.nStates];
		betaNorms = new double[modelInfo.longestSeq];
	}
	
	public double apply(double[] param, double[] grad) {
		// Initialize values
		Arrays.fill(grad, 0);
		double result = 0.0;
		Arrays.fill(expects, 0);

		// Iterate through sequences
		int seqStart = 0;
		for(int i = 0; i<modelInfo.nSeqs; ++i) {
			int len = modelInfo.seqOffsets[i+1] - modelInfo.seqOffsets[i]; 
			
			// Work backwards, computing betas.
			Arrays.fill(betas[len-1], 1.0);
			betaNorms[len-1] = 0;
			for (int pos = len - 1; pos > 0; --pos) {
				calcMi(i, pos, param);
				quickBetaUpdate(betas[pos], betas[pos-1]);

				double n = normalizePotential(betas[pos-1]);
				betaNorms[pos-1] = betaNorms[pos] + log(n);
			}
			
			// Now work forwards
			double logZ = Double.NEGATIVE_INFINITY;	// This should always get initialized.  Blow up if not.
			double alphaNorm = 0;
			double prevAlphaNorm = 0;
			for(int pos=0; pos<len; ++pos) {
				double[] beta = betas[pos];
				double betaNorm = betaNorms[pos];
				if(pos == 0) {
					calcStartAlpha(i, param);
					alphaNorm = log(normalizePotential(alpha));
					
					// We now have everything needed to compute Z.
					logZ = log(ColtUtil.dotProduct(alpha, beta)) + betaNorm + alphaNorm;
					//log.info("Z = "+logZ);
				}
				else {
					calcMi(i, pos, param);
					//r * M (or M'*r)
					//mi.zMult(prevAlpha, alpha, 1, 0, true);
					quickAlphaUpdate(prevAlpha, alpha);
					alphaNorm = prevAlphaNorm + log(normalizePotential(alpha));
					
					// Verify our calculations by checking the Z.
					// Expensive, so only enable as needed during debugging
					double newZ = log(ColtUtil.dotProduct(alpha, beta)) + betaNorm + alphaNorm;
					Assert.a(Math.abs(newZ-logZ) < 0.0000001*Math.abs(logZ), "New Z:",newZ," Old was: ", logZ);
				}
				
				// Iterate through the potentials to update feature expectations
				double nodeNorm = exp(alphaNorm + betaNorm - logZ);
				double edgeNorm = exp(prevAlphaNorm + betaNorm - logZ);
				updateExpectations(i, pos, nodeNorm, edgeNorm, beta);

				if(debug) {
					if((i < 2 || i == modelInfo.nSeqs-1) && (pos < 2 || pos >= len-2)) {
						log.debug(String.format("Pos: %d expects: %s alphas: %s (norm %f) betas: %s (norm %f)", pos, ColtUtil.format(expects), ColtUtil.format(alpha), alphaNorm, ColtUtil.format(beta), betaNorm));
					}
				}

				// Recycle the arrays
				double[] swap = prevAlpha;
				prevAlpha = alpha;
				alpha = swap;
				prevAlphaNorm = alphaNorm;
			}

			result -=  logZ;
			seqStart += len;
		}
		
		// sum_j lambda_j F_j(xk, yk)
		double[] featureSums = cacheProcessor.getFeatureSums();
		for (int j = 0; j < modelInfo.nFeatures; ++j) {
			result += featureSums[j] * param[j];
			grad[j] = featureSums[j] - expects[j];
		}

		if(log.isInfoEnabled()) {
			// Report average per-label numbers.
			log.info(String.format("It: %d L=%e, LL=%f, norm(grad): %f Sums: %s Expects: %s Weights: %s Grad (unnorm): %s", iter, exp(result/modelInfo.totalPositions), result/modelInfo.totalPositions, ColtUtil.norm(grad)/modelInfo.totalPositions, ColtUtil.format(featureSums), ColtUtil.format(expects), ColtUtil.format(param), ColtUtil.format(grad)));
		}
		iter += 1;
		result = result/modelInfo.totalPositions; 
		for(int i=0; i<grad.length; ++i) {
			grad[i] = grad[i]/modelInfo.totalPositions;
		}
	
		this.featureSums = featureSums;

		return result;
	}

	public void clean() {
	}
	
	/** This is one of the most time critical parts of the entire solver.  The goal is to update the transition matrix.
	 * This function makes a lot of assumptions in order to maximize performance.
	 * 
	 * To maximize performance, we want to make one pass through the Mi matrix, setting each entry to its correct value.
	 * The value for each entry is the exponent of the sum of the weighted feature values of the edge for that entry and 
	 * its corresponding node.  The entry s0,s1 consists of the s0,s1 edge and the s1 node. 
	 * 
	 * Because node features are applied to more than 1 entry in the matrix, we use a sorting of all of the features where
	 * each node preceeds all its corresponding edges.  This allows us to keep track of only 1 node value at a time and easily
	 * apply it to all its edge features.
	 * 
	 * As we evaluate each potential, we check the cache to see if it is valid at this position and to get any features values.
	 * Note that this function very much depends on the fact that the entries in the cache will be in the correct order.
	 * 
	 * The other wrinkle is that for features that always occur (constant features), we pull them from the constant mi array, not
	 * from the cache.
	 * 
	 *  This function is also used to calculate the Mi Matrix for the constant features.
	 *  */
	void calcMi(int seq, int pos, double[] lambda) {
		cacheProcessor.evaluatePosition(seq, pos);
		double nodeVal = Double.NaN;
		int overallPosition = modelInfo.seqOffsets[seq]+pos;
		int invalidIndex = overallPosition*modelInfo.nPotentials;
		for(short potential : modelInfo.orderedPotentials) {
			boolean invalid = invalidTransitions[invalidIndex + potential];
			double features = 0.0;
			
			if(invalid) {
				features = Double.NEGATIVE_INFINITY;
			}
			else {
				// Add up all features for this potential.
				FeatureEvaluation potEvals = evals[potential];
				short[] indices = potEvals.index;
				float[] vals = potEvals.value;
				int i = 0;
				short index = indices[i];
				while(index >= 0) {
					features += vals[i]*lambda[index]; 
					index = indices[++i];
				}
				if(index == Short.MIN_VALUE) {
					// An invalid potential is indicated by a feature value of Short.MIN_VALUE
					features = Double.NEGATIVE_INFINITY; 
				}
			}
			
			if(potential < modelInfo.nStates) {
				nodeVal = features;
			}
			else {
				int transition = potential - modelInfo.nStates;
				mi[transition] = exp(features + nodeVal);
			}
		}		
	}
	
	/** A specialized version of calcMi for the first position in a sequence. Has the special property that 
	 * constant edge features are not included.  Also optionally allows you to set a node value used to initialize the alphas. */
	void calcStartAlpha(int seq, double[] lambda) {
		cacheProcessor.evaluatePosition(seq, 0);
		int overallPosition = modelInfo.seqOffsets[seq];
		int invalidIndex = overallPosition*modelInfo.nPotentials;
		for(short potential : modelInfo.orderedPotentials) {
			if(potential < modelInfo.nStates) {
				boolean invalid = invalidTransitions[invalidIndex + potential];
				double features = 0.0;

				if(invalid) {
					features = Double.NEGATIVE_INFINITY;
				}
				else {
					// Add up all features for this potential.
					FeatureEvaluation potEvals = evals[potential];
					short[] indices = potEvals.index;
					float[] vals = potEvals.value;
					int i = 0;
					short index = indices[i];
					while(index >= 0) {
						features += vals[i]*lambda[index];
						index = indices[++i];
					}
					if(index == Short.MIN_VALUE) {
						features = Double.NEGATIVE_INFINITY;
					}
				}
				alpha[potential] = exp(features);
			}
		}		
	}
	
	private void quickBetaUpdate(double[] lastBeta, double[] newBeta) {
		Arrays.fill(newBeta, 0);
		double nodeVal = 0.0;
		for(short potential : modelInfo.orderedPotentials) {
			if(potential < modelInfo.nStates) {
				nodeVal = lastBeta[potential];
			}
			else {
				int trans = potential - modelInfo.nStates;
				int from = modelInfo.transitionFrom[trans];
				newBeta[from] += mi[trans]*nodeVal;
				//log.debug(String.format("beta[%d] = %f = %f + mi: %f * last: %f",
				//			from, newBeta[from], newBeta[from]-mi[trans]*nodeVal, mi[trans], nodeVal));
			}
		}
	}

	private void quickAlphaUpdate(double[] lastAlpha, double[] newAlpha) {
		double nodeVal = 0.0;
		int lastState = -1;
		for(short potential : modelInfo.orderedPotentials) {
			if(potential < modelInfo.nStates) {
				if(lastState != -1) {
					newAlpha[lastState] = nodeVal;
				}
				lastState = potential;
				nodeVal = 0.0;
			}
			else {
				int trans = potential - modelInfo.nStates;
				int from = modelInfo.transitionFrom[trans];
				nodeVal += lastAlpha[from]*mi[trans];
			}
		}
		newAlpha[lastState] = nodeVal;
	}

	void updateExpectations(int seq, int pos, double nodeNorm, double edgeNorm, double[] beta) {
		int currentNode = -1;
		double currentBeta = 0.0f;
		int overallPos = modelInfo.seqOffsets[seq]+pos;
		int invalidIndex = overallPos*modelInfo.nPotentials;
		for(short potential : modelInfo.orderedPotentials) {
			boolean invalid = invalidTransitions[invalidIndex + potential];
			if(invalid)
				continue;
			
			double prob = 0.0;
			if(potential < modelInfo.nStates) {
				currentNode = potential;
				currentBeta = beta[currentNode];
				//log.info(String.format("State prob #%d: %f * %f * %f ", currentNode, alpha[currentNode], currentBeta, nodeNorm));
				prob = alpha[currentNode] * currentBeta * nodeNorm;
			}
			else {
				if(pos == 0)
					continue;
				int trans = potential - modelInfo.nStates;
				int yprev = modelInfo.transitionFrom[trans];
				prob = prevAlpha[yprev] * mi[trans] * currentBeta * edgeNorm;
			}

			// Iterate through features for this potential.
			FeatureEvaluation potEvals = evals[potential];
			short[] indices = potEvals.index;
			float[] vals = potEvals.value;
			int i = 0;
			short index = indices[i];
			if(index != Short.MIN_VALUE) {
				while(index != -1) {
					//log.info(String.format("Expect #%d: %f * %f ", index, prob, vals[i]));
					expects[index] += prob*vals[i];
					index = indices[++i];
				}
			}
		}
	}
	
	private double normalizePotential(double[] vec) {
		double norm = 0.0;
		int len = vec.length;
		for(int i=0; i<len; ++i) {
			norm += vec[i];
		}
		double mult = 1/norm;
		for(int i=0; i<len; ++i) {
			vec[i] *= mult;
		}
		return norm;
	}

	static final double exp(double val) {
		return Math.exp(val);
	}
	
	static final double log(double val) {
		return Math.log(val);
	}

	public double[] getFeatureSums() {
		return this.featureSums.clone();
	}
}