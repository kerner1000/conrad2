/**
 * 
 */
package calhoun.analysis.crf.solver.check;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.CRFObjectiveFunctionGradient;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.analysis.crf.solver.CacheProcessor;
import calhoun.util.Assert;
import calhoun.util.ColtUtil;

public class OldCachedCRFGradient implements CRFObjectiveFunctionGradient {
	private static final Log log = LogFactory.getLog(OldCachedCRFGradient.class);
	boolean debug = log.isDebugEnabled();

	boolean allPaths;
	
	/// Cache feature information
	short[] id;
	byte[] potentialIx;
	float[] val;

	/// Cached value of the Mi matrix for all of the features present at every position
	/// Mi is stored as a sparse matrix
	int miLength;
	double[] constMi;
	double[] mi;
	short[] transitionFrom;
	short[] transitionTo;
	short[] orderedPotentials;
	boolean[] invalidTransitions;
	int totalPositions;

	/// Cached values of the sums of each feature value through the whole training set.
	double[] featureSums;

	/// Index into the feature arrays of the first feature for each postion of each sequence.
	int[] starts;

	/// Index into the starts array of the first position of each sequence.
	int[] seqOffsets;

	/// Number of sequences in the training data set
	int nSeqs; 
	
	/// Number of constant features.
	int nConstantFeatures; 
	
	List<? extends TrainingSequence<?>> data;
	ModelManager fm;
	int nFeatures;
	int nStates;
	int nPotentials;
	int iter = 0;

	double[] prevAlpha;
	double[] alpha;
	double[][] betas;
	double[] betaNorms;
	double[] expects;

	double exp(double val1) {
		return Math.exp(val1);
	}
	
	double log(double val1) {
		return Math.log(val1);
	}
	
	public OldCachedCRFGradient(boolean allPaths) {
		this.allPaths = allPaths;
	}

	public void setTrainingData(ModelManager fm, List<? extends TrainingSequence<?>> data) { 
		this.fm = fm;
		this.data = data;
		nFeatures = fm.getNumFeatures();
		nStates = fm.getNumStates();
		nSeqs = data.size();

		expects = new double[nFeatures];
		prevAlpha = new double[nStates];
		alpha = new double[nStates];

		FeatureCache cache = new FeatureCache(fm, data, allPaths);

		// Create local references to cache values
		orderedPotentials = cache.orderedPotentials;
		id = cache.id;
		potentialIx = cache.potentialIx;
		val = cache.val;

		transitionFrom = cache.transitionFrom;
		transitionTo = cache.transitionTo;
		miLength = cache.nTransitions;
		constMi = new double[miLength];
		mi = new double[miLength];

		featureSums = cache.featureSums;
		starts = cache.starts;
		seqOffsets = cache.seqOffsets;
		invalidTransitions = cache.invalidTransitions;
		nConstantFeatures = cache.numConstantFeatures;
		nPotentials = cache.nPotentials;
		totalPositions = cache.totalPositions;
		
		// Initialize betas (for use later, in the gradient computation)
		betas = new double[cache.longestSeq][nStates];
		betaNorms = new double[cache.longestSeq];
	}
	
	public void clean() {
	}
	
	public double apply(double[] param, double[] grad) {
		Arrays.fill(grad, 0);
		double result = 0.0;

		// Calculate the constant Mi matrix
		// CalcMi produces an exponentiated matrix, so we take the log the get the right init values
		// and copy them to the constMi matrix.
		Arrays.fill(constMi, 0.0);
		calcMi(-1, 0, starts[0], param);
		for(int i = 0; i<miLength; ++i) {
			constMi[i] = log(mi[i]);
		}
		
		// Iterate through sequences
		Arrays.fill(expects, 0);
		int seqStart = 0;
		for(int i = 0; i<nSeqs; ++i) {
			int len = seqOffsets[i+1] - seqOffsets[i]; 
			
			// Work backwards, computing betas.
			Arrays.fill(betas[len-1], 1.0);
			betaNorms[len-1] = 0;
			int cacheStop = starts[seqStart + len];
			for (int pos = len - 1; pos > 0; --pos) {
				int overallPosition = seqStart + pos;
				int cacheStart = starts[overallPosition];
				calcMi(overallPosition, cacheStart, cacheStop, param);
				cacheStop = cacheStart;
				quickBetaUpdate(betas[pos], betas[pos-1]);
				
				//mi.zMult(betas[j], betas[j-1], 1, 0, false);
				double n = normalizePotential(betas[pos-1]);
				betaNorms[pos-1] = betaNorms[pos] + log(n);
			}
			
			// Now work forwards
			double logZ = Double.NEGATIVE_INFINITY;	// This should always get initialized.  Blow up if not.
			double alphaNorm = 0;
			double prevAlphaNorm = 0;
			int cacheStart = starts[seqStart];
			for(int pos=0; pos<len; ++pos) {
				int overallPosition = seqStart + pos;
				double[] beta = betas[pos];
				double betaNorm = betaNorms[pos];
				cacheStop = starts[overallPosition+1];
				if(pos == 0) {
					calcStartAlpha(overallPosition, cacheStart, cacheStop, param);
					alphaNorm = log(normalizePotential(alpha));
					
					// We now have everything needed to compute Z.
					logZ = log(ColtUtil.dotProduct(alpha, beta)) + betaNorm + alphaNorm;
					//log.info("Z = "+logZ);
				}
				else {
					calcMi(overallPosition, cacheStart, cacheStop, param);
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
				updateExpectations(overallPosition, pos != 0, cacheStart, cacheStop, nodeNorm, edgeNorm, beta);

				if(debug) {
					if((i < 2 || i == nSeqs-1) && (pos < 2 || pos >= len-2)) {
						log.debug(String.format("Pos: %d expects: %s alphas: %s (norm %f) betas: %s (norm %f)", pos, ColtUtil.format(expects), ColtUtil.format(alpha), alphaNorm, ColtUtil.format(beta), betaNorm));
					}
				}

				// Recycle the arrays
				double[] swap = prevAlpha;
				prevAlpha = alpha;
				alpha = swap;
				prevAlphaNorm = alphaNorm;
				cacheStart = cacheStop;
			}

			result -=  logZ;
			seqStart += len;
		}
		
		// sum_j lambda_j F_j(xk, yk)
		for (int j = 0; j < nFeatures; ++j) {
			result += featureSums[j] * param[j];
			grad[j] = featureSums[j] - expects[j];
		}

		if(log.isInfoEnabled()) {
			// Report average per-label numbers.
			log.info(String.format("It: %d L=%e, LL=%f, norm(grad): %f Sums: %s Expects: %s Weights: %s Grad (unnorm): %s", iter, exp(result/totalPositions), result/totalPositions, ColtUtil.norm(grad)/totalPositions, ColtUtil.format(featureSums), ColtUtil.format(expects), ColtUtil.format(param), ColtUtil.format(grad)));
			//log.info(String.format("It: %d L=%e, LL=%f, norm(grad): %f Weights: %s Grad: %s", iter, exp(result), result, ColtUtil.norm(grad), ColtUtil.format(param), ColtUtil.format(grad)));
		}
		iter += 1;
		result = result/totalPositions; 
		for(int i=0; i<grad.length; ++i) {
			grad[i] = grad[i]/totalPositions;
		}
		return result;
	}

	void updateExpectations(int overallPos, boolean includeEdges, int posCurrent, int posStop, double nodeNorm, double edgeNorm, double[] beta) {
		int constCurrent = 0;
		
		// Constant features
		short constId = -1;
		byte constPotential = -1;
		double constVal = Double.NaN;
		if(constCurrent < nConstantFeatures) {
			constId = id[constCurrent];
			constPotential = potentialIx[constCurrent];
			constVal = val[constCurrent];
			++constCurrent;
		}
		
		// Positional features
		short posId = -1;
		byte posPotential = -1;
		double posVal = Double.NaN;
		if(posCurrent < posStop) {
			posId = id[posCurrent];
			posPotential = potentialIx[posCurrent];
			posVal = val[posCurrent];
			++posCurrent;
		}

		int currentNode = -1;
		double currentBeta = 0.0f;
		int invalidIndex = overallPos*nPotentials;
		for(short potential : orderedPotentials) {
			boolean invalid = invalidTransitions[invalidIndex + potential];

			double prob = 0.0;
			if(potential < nStates) {
				currentNode = potential;
				currentBeta = beta[currentNode];
				if(!invalid) {
					//log.info(String.format("State prob #%d: %f * %f * %f ", currentNode, alpha[currentNode], currentBeta, nodeNorm));
					prob = alpha[currentNode] * currentBeta * nodeNorm;
				}
			}
			else {
				int trans = potential - nStates;
				int yprev = transitionFrom[trans];
				if(!invalid) {
					prob = prevAlpha[yprev] * mi[trans] * currentBeta * edgeNorm;
				}
			}
			
			// Include constant features for this potential.  Need to skip over constant features for invalid potentials
			while(constPotential == potential) {
				if(!invalid && (includeEdges || potential < nStates)) {
					//log.info(String.format("Expect #%d: %f * %f ", constId, prob, constVal));
					expects[constId] += prob*constVal;
				}
				if(constCurrent < nConstantFeatures) {
					constId = id[constCurrent];
					constVal = val[constCurrent];
					constPotential = potentialIx[constCurrent];
					++constCurrent;
				}
				else { break; }
			}

			// Include cached features for this potential.  The cache shoudl never have features for invalid potentials.
			if(!invalid) {
				while(posPotential == potential) {
					//log.info(String.format("Expect #%d: %f * %f ", posId, prob, posVal));
					expects[posId] += prob*posVal;
					if(posCurrent < posStop) {
						posId = id[posCurrent];
						posVal = val[posCurrent];
						posPotential = potentialIx[posCurrent];
						++posCurrent;
					}
					else { break; }
				}
			}
		}
		// verify that we actually checked all the cache entries.
		Assert.a(constCurrent == nConstantFeatures);
		Assert.a(posCurrent == posStop);
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
	void calcMi(int overallPosition, int current, int stop, double[] lambda) {
		// Features will always be in the order of the potentials.  We loop through the potentials, grabbing the features for each.
		byte cachedPotential = -1;
		double cachedVal = Double.NaN;
		if(current < stop) {
			cachedPotential = potentialIx[current];
			cachedVal = val[current]*lambda[id[current]];
			++current;
		}
		double nodeVal = Double.NaN;
		int invalidIndex = overallPosition*nPotentials;
		for(short potential : orderedPotentials) {
			boolean invalid = overallPosition != -1 && invalidTransitions[invalidIndex + potential];
			double features = invalid ? Double.NEGATIVE_INFINITY : 0.0;

			// Add up all features for this potential.
			while(cachedPotential == potential) {
				features += cachedVal;
				//Assert.a(!Double.isNaN(features));
				if(current < stop) {
					cachedVal = val[current]*lambda[id[current]];
					cachedPotential = potentialIx[current];
					++current;
				}
				else { break; }
			}

			if(potential < nStates) {
				nodeVal = features;
			}
			else {
				//log.debug(String.format("Mi[%d, %d] = %f, adding in %f to get %f", feat.yprev(), feat.y(), val, feat.value(), val*exp(feat.value()*param[feat.index()])));
				int transition = potential - nStates;
				mi[transition] = exp(features + nodeVal + constMi[transition]);
			}
		}
		// verify that we actually checked all the cache entries.
		Assert.a(current == stop, "Pos: ", overallPosition, " Expected ", stop, " features only found ", current);
	}
	
	/** A specialized version of calcMi for the first position in a sequence. Has the special property that 
	 * constant edge features are not included.  Also optionally allows you to set a node value used to initialize the alphas. */
	void calcStartAlpha(int overallPosition, int posCurrent, int posStop, double[] lambda) {
		int constCurrent = 0;
		
		// Constant features
		byte constPotential = -1;
		double constVal = Double.NaN;
		if(constCurrent < nConstantFeatures) {
			constPotential = potentialIx[constCurrent];
			constVal = val[constCurrent]*lambda[id[constCurrent]];
			++constCurrent;
		}
		
		// Positional features
		byte posPotential = -1;
		double posVal = Double.NaN;
		if(posCurrent < posStop) {
			posPotential = potentialIx[posCurrent];
			posVal = val[posCurrent]*lambda[id[posCurrent]];
			++posCurrent;
		}

		int invalidIndex = overallPosition*nPotentials;
		for(short potential : orderedPotentials) {
			boolean invalid = invalidTransitions[invalidIndex + potential];
			double features = invalid ? Double.NEGATIVE_INFINITY : 0.0;

			// Find constant features for this potential.  Ignore edge features.
			while(constPotential == potential) {
				if(potential < nStates) {
					features += constVal;
					Assert.a(!Double.isNaN(features));
				}
				if(constCurrent < nConstantFeatures) {
					constVal = val[constCurrent]*lambda[id[constCurrent]];
					constPotential = potentialIx[constCurrent];
					++constCurrent;
				}
				else { break; }
			}

			if(potential < nStates) {
				// Add in cached features
				while(posPotential == potential) {
					features += posVal;
					if(posCurrent < posStop) {
						posVal = val[posCurrent]*lambda[id[posCurrent]];
						posPotential = potentialIx[posCurrent];
						++posCurrent;
					}
					else { break; }
				}
				alpha[potential] = exp(features);
			}
		}
		// verify that we actually checked all the cache entries.
		Assert.a(constCurrent == nConstantFeatures);
		Assert.a(posCurrent == posStop);
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

	private void quickBetaUpdate(double[] lastBeta, double[] newBeta) {
		Arrays.fill(newBeta, 0);
		double nodeVal = 0.0;
		for(short potential : orderedPotentials) {
			if(potential < nStates) {
				nodeVal = lastBeta[potential];
			}
			else {
				int trans = potential - nStates;
				int from = transitionFrom[trans];
				newBeta[from] += mi[trans]*nodeVal;
			}
		}
	}

	private void quickAlphaUpdate(double[] lastAlpha, double[] newAlpha) {
		double nodeVal = 0.0;
		int lastState = -1;
		for(short potential : orderedPotentials) {
			if(potential < nStates) {
				if(lastState != -1) {
					newAlpha[lastState] = nodeVal;
				}
				lastState = potential;
				nodeVal = 0.0;
			}
			else {
				int trans = potential - nStates;
				int from = transitionFrom[trans];
				nodeVal += lastAlpha[from]*mi[trans];
			}
		}
		newAlpha[lastState] = nodeVal;
	}

	public void setCacheProcessor(CacheProcessor cacheProcessor) {
		throw new UnsupportedOperationException("This objective function predates cacheProcessors");
	}
}