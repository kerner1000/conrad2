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
import calhoun.util.ColtUtil;
import calhoun.util.Util;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Blas;
import cern.colt.matrix.linalg.SeqBlas;

public class BasicCRFGradient implements CRFObjectiveFunctionGradient {
	private static final Log log = LogFactory.getLog(BasicCRFGradient.class);
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
	DoubleMatrix1D expects;
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
	}

	public void clean() {
	}
	
	public double apply(double[] param, double[] grad) {
		TransitionInfo t = new TransitionInfo(fm, allPaths);
		FeatureCalculator calc = new FeatureCalculator(fm, param, t);
		Arrays.fill(grad, 0);
		int seqIx = 0;
		double result = 0.0;
		// Allocate the mi and feature sum matrices
		for (TrainingSequence seq : data) {
			int len = seq.length();
			// Initialize
			expects.assign(0);

			// Resize the beta array so it fits the largest sequence
			if(betas.length < len) {
				betas = new DenseDoubleMatrix1D[len];
				for(int i = 0; i<betas.length; ++i) {
					betas[i] = new DenseDoubleMatrix1D(nStates);
				}
			}

			// Work backwards, computing betas.
			betas[len-1].assign(1.0);
			for (int i = len - 1; i > 0; --i) {
				// alpha is just used as a temp here
				mi.assign(0);
				calc.computeMi(seq, i, mi, null);
				mi.assign(ColtUtil.exp);
				mi.zMult(betas[i], betas[i-1], 1, 0, false);
			}

			calc.resetFeatureSums();
			
			// Now work forwards
			for(int pos=0; pos<len; ++pos) {
				if(pos == 0) {
					for(int state = 0; state < nStates; ++state) {
						double val = calc.calcNodeValue(seq, pos, state);
						alpha.setQuick(state, val);
					}
					alpha.assign(ColtUtil.exp);
				}
				else {
					//r * M (or M'*r)
					mi.assign(0);
					calc.computeMi(seq, pos, mi, alpha);
					mi.assign(ColtUtil.exp);
					mi.zMult(prevAlpha, alpha, 1, 0, true);
				}

				// Iterate through the features to update the expectation for each one
				ArrayFeatureList results = new ArrayFeatureList(fm);
				for(int state = 0; state < nStates; ++state) {
					results.evaluateNode(seq, pos, state);
					double mult = alpha.getQuick(state) * betas[pos].getQuick(state);
					results.updateExpectations(expects, mult);

					if(pos > 0) {
						for(int prevState = 0; prevState < nStates; ++prevState) {
							if (calc.isValidTransition(prevState, state)) {
								results.evaluateEdge(seq, pos, prevState, state);
								mult = prevAlpha.getQuick(prevState) * mi.getQuick(prevState, state) * betas[pos].getQuick(state);
								results.updateExpectations(expects, mult);
							}
						}
					}
				}

				// Recycle the arrays
				DoubleMatrix1D swap = prevAlpha;
				prevAlpha = alpha;
				alpha = swap;

				//log.debug(String.format("Pos: %d expects: %s alphas: %s betas: %s", pos, ColtUtil.format(expects.toArray()), ColtUtil.format(alpha.toArray()), ColtUtil.format(betas[pos].toArray())));
			}

			// sum_j lambda_j F_j(xk, yk)
			double z = prevAlpha.zSum();
			double seqLogPotential = calc.getWeightedFeatureSum();
			double[] featureSums = calc.getFeatureSums();
			for (int i = 0; i < nFeatures; ++i) {
				grad[i] += featureSums[i] - expects.getQuick(i)/z;
			}
			result += seqLogPotential - Math.log(z);
			if(debug) {
				log.debug(String.format("Seq: %d L=%e, LL=%f, Feats=%s, Sum=%e (log=%f), Z=%e (log=%f) Expects=%s, Grad=%s", seqIx, Math.exp(seqLogPotential - Math.log(z)), seqLogPotential - Math.log(z),
					StringUtils.join(Util.convertDoubleArray(featureSums).iterator(), ','), Math.exp(seqLogPotential), seqLogPotential, z,
					Math.log(z), ColtUtil.format(expects.toArray()), ColtUtil.format(grad)));
			}
			seqIx += 1;
		}
		if(log.isInfoEnabled()) {
			log.info(String.format("It: %d L=%e, LL=%f, norm(grad): %f Weights: %s Grad: %s", iter, Math.exp(result), result, ColtUtil.norm(grad), ColtUtil.format(param), ColtUtil.format(grad)));
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

	public boolean isAllPaths() {
		return allPaths;
	}

	public void setAllPaths(boolean allPaths) {
		this.allPaths = allPaths;
	}
}