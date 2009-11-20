/**
 * 
 */
package calhoun.analysis.crf.solver.check;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.lang.StringUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.CRFObjectiveFunctionGradient;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.analysis.crf.solver.CacheProcessor;
import calhoun.util.Assert;
import calhoun.util.ColtUtil;
import calhoun.util.Util;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Blas;
import cern.colt.matrix.linalg.SeqBlas;
import cern.jet.math.Mult;

public class NormalizedCRFGradient implements CRFObjectiveFunctionGradient {
	private static final Log log = LogFactory.getLog(NormalizedCRFGradient.class);
	boolean debug = log.isDebugEnabled();
	
	List<? extends TrainingSequence<?>> data;
	ModelManager fm;
	int nFeatures;
	int nStates;
	double[] beta = null;
	Blas blas = SeqBlas.seqBlas;
	int iter = 0;
	DoubleMatrix2D mi;
	DoubleMatrix1D ri;
	DoubleMatrix1D temp;
	DoubleMatrix1D prevAlpha;
	DoubleMatrix1D alpha;
	DoubleMatrix1D[] betas;
	double[] betaNorm;
	DoubleMatrix1D expects;
	Mult normalizer = Mult.mult(0);
	boolean allPaths = false;
	
	public void setTrainingData(ModelManager fm, List<? extends TrainingSequence<?>> data) { 
		this.fm = fm;
		this.data = data;
		nFeatures = fm.getNumFeatures();
		nStates = fm.getNumStates();
		expects = new DenseDoubleMatrix1D(nFeatures);
		mi = new DenseDoubleMatrix2D(nStates, nStates);
		prevAlpha = new DenseDoubleMatrix1D(nStates);
		alpha = new DenseDoubleMatrix1D(nStates);
		betas = new DenseDoubleMatrix1D[0];
		betaNorm = new double[0];
	}

	double normalizePotential(DoubleMatrix1D vec) {
		double norm = vec.zSum();
		normalizer.multiplicator = 1/norm;
		vec.assign(normalizer);
		return norm;
	}
	
	public void clean() {
	}
	
	public double apply(double[] param, double[] grad) {
		TransitionInfo t = new TransitionInfo(fm, allPaths);
		FeatureCalculator calc = new FeatureCalculator(fm, param, t);
		Arrays.fill(grad, 0);
		int seqIx = 0;
		double result = 0.0;
		double[] totalFeatureSums = new double[nFeatures];
		double[] totalExpects = new double[nFeatures];
		// Allocate the mi and feature sum matrices
		for (TrainingSequence seq : data) {
			int len = seq.length();

			// Resize the beta array so it fits the largest sequence
			if(betas.length < len) {
				betas = new DenseDoubleMatrix1D[len];
				betaNorm = new double[len];
				for(int i = 0; i<betas.length; ++i) {
					betas[i] = new DenseDoubleMatrix1D(nStates);
				}
			}

			// Work backwards, computing betas.
			// TODO: betas for nodes that are invalid at the last position should be reset to 0. 
			betas[len-1].assign(1.0);
			betaNorm[len-1] = 0;
			for (int i = len - 1; i > 0; --i) {
				mi.assign(0);
				calc.computeMi(seq, i, mi, null);
				mi.assign(ColtUtil.exp);
				mi.zMult(betas[i], betas[i-1], 1, 0, false);
				double n = normalizePotential(betas[i-1]);
				betaNorm[i-1] = betaNorm[i] + Math.log(n);
			}
			
			// Now work forwards
			// Initialize
			calc.resetFeatureSums();
			expects.assign(0);

			double logZ = Double.NEGATIVE_INFINITY;	// This should always get initialized.  Blow up if not.
			double alphaNorm = 0;
			double prevAlphaNorm = 0;
			for(int pos=0; pos<len; ++pos) {
				if(pos == 0) {
					for(int state = 0; state < nStates; ++state) {
						double val = calc.calcNodeValue(seq, pos, state);
						alpha.setQuick(state, val);
					}
					alpha.assign(ColtUtil.exp);
					alphaNorm = Math.log(normalizePotential(alpha));
					
					// We now have everything needed to compute Z.
					logZ = Math.log(alpha.zDotProduct(betas[0])) + betaNorm[0] + alphaNorm;
				}
				else {
					mi.assign(0);
					calc.computeMi(seq, pos, mi, alpha);
					mi.assign(ColtUtil.exp);

					//r * M (or M'*r)
					mi.zMult(prevAlpha, alpha, 1, 0, true);
					alphaNorm = prevAlphaNorm + Math.log(normalizePotential(alpha));

					// Verify our calculations by checking the Z.
					double newZ = Math.log(alpha.zDotProduct(betas[pos])) + betaNorm[pos] + alphaNorm;
					Assert.a(Math.abs(newZ-logZ) < 0.0000001*Math.abs(logZ), "New Z:",newZ," Old was: ", logZ);
				}

				// Iterate through the features to update the expectation for each one
				ArrayFeatureList results = new ArrayFeatureList(fm);
				double nodeNorm = Math.exp(alphaNorm + betaNorm[pos] - logZ);
				double edgeNorm = Math.exp(prevAlphaNorm + betaNorm[pos] - logZ);
				for(int state = 0; state < nStates; ++state) {
					results.evaluateNode(seq, pos, state);
					double mult = alpha.getQuick(state) * betas[pos].getQuick(state) * nodeNorm;
					results.updateExpectations(expects, mult);

					if(pos > 0) {
						for(int prevState = 0; prevState < nStates; ++prevState) {
							if (calc.isValidTransition(prevState, state)) {
								results.evaluateEdge(seq, pos, prevState, state);
								mult = prevAlpha.getQuick(prevState) * mi.getQuick(prevState, state) * betas[pos].getQuick(state) * edgeNorm;
								results.updateExpectations(expects, mult);
							}
						}
					}
				}

				if(debug) {
					if((seqIx < 2 || seqIx == data.size()-1) && (pos < 2 || pos >= len-2)) {
						log.debug(String.format("Pos: %d expects: %s alphas: %s betas: %s", pos, ColtUtil.format(expects.toArray()), ColtUtil.format(alpha.toArray()), ColtUtil.format(betas[pos].toArray())));
					}
				}

				// Recycle the arrays
				DoubleMatrix1D swap = prevAlpha;
				prevAlpha = alpha;
				alpha = swap;
				prevAlphaNorm = alphaNorm;
			}

			// sum_j lambda_j F_j(xk, yk)
			double seqLogLikelihood = calc.getWeightedFeatureSum();
			double[] featureSums = calc.getFeatureSums();
			for (int i = 0; i < nFeatures; ++i) {
				grad[i] += featureSums[i] - expects.getQuick(i);
				totalFeatureSums[i] += featureSums[i];
				totalExpects[i] += expects.getQuick(i);
			}
			result += seqLogLikelihood - logZ;
			if(debug) {
				log.debug(String.format("Seq: %d L=%e, LL=%f, Feats=%s, Sum=%e (log=%f), Z=%e (log=%f) Expects=%s, Grad=%s", seqIx, Math.exp(seqLogLikelihood - logZ), seqLogLikelihood - logZ,
					StringUtils.join(Util.convertDoubleArray(featureSums).iterator(), ','), Math.exp(seqLogLikelihood), seqLogLikelihood, Math.exp(logZ),
					logZ, ColtUtil.format(expects.toArray()), ColtUtil.format(grad)));
			}
			seqIx += 1;
		}

		// Make the result average likelihood per label
		if(log.isInfoEnabled()) {
			log.info(String.format("It: %d L=%e, LL=%f, norm(grad): %f Sums: %s Expects: %s Weights: %s Grad: %s", iter, Math.exp(result), result, ColtUtil.norm(grad), ColtUtil.format(totalFeatureSums), ColtUtil.format(totalExpects), ColtUtil.format(param), ColtUtil.format(grad)));
			//log.info(String.format("It: %d L=%e, LL=%f, norm(grad): %f Weights: %s Grad: %s", iter, Math.exp(result), result, ColtUtil.norm(grad), ColtUtil.format(param), ColtUtil.format(grad)));
		}
		iter += 1;
		int totalPositions = 0;
		for(TrainingSequence i : data) {
			totalPositions += i.length();
		}
		result = result/totalPositions; 
		for(int i=0; i<grad.length; ++i) {
			grad[i] = grad[i]/totalPositions;
		}
		return result;
	}

	public void setCacheProcessor(CacheProcessor cacheProcessor) {
		throw new UnsupportedOperationException("This objective function predates cacheProcessors");
	}

	public boolean isAllPaths() {
		return allPaths;
	}

	public void setAllPaths(boolean allPaths) {
		this.allPaths = allPaths;
	}
}