package calhoun.analysis.crf.solver.semimarkov;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.CRFObjectiveFunctionGradient;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.analysis.crf.solver.CacheProcessor;
import calhoun.analysis.crf.solver.LogFiles;
import calhoun.analysis.crf.solver.LookbackBuffer;
import calhoun.analysis.crf.solver.RecyclingBuffer;
import calhoun.analysis.crf.solver.CacheProcessor.FeatureEvaluation;
import calhoun.analysis.crf.solver.CacheProcessor.LengthFeatureEvaluation;
import calhoun.analysis.crf.solver.CacheProcessor.SolverSetup;
import calhoun.analysis.crf.solver.CacheProcessor.StatePotentials;
import calhoun.util.Assert;
import calhoun.util.ColtUtil;

/** computes the likelihood of the true path for a semi-Markov CRF.  The likelihood is normalized to a per label likelihood. 
 * <h2>Debugging output</h2>
 * To get a better understanding of what the objective function is doing, several differn properties can be set that
 * cause the objective function to write out trace files showing its calculations during training.  Usually when turning
 * these options on, you should set <code>maxIters = 1</code> and <code>requireConvergence = false</code> in your optimizer
 * to do only a single training iteration, possibly setting the starts to some predetermined value.  Each of these
 * properties can be configured with a filename and each time {@link #apply} is called, the file will be overwritten with 
 * data from the current call.  The logging options are:
 * <ul>
 * <li> <b><code>alphaFile</code></b> - computation of alpha values for Markov states, includes all nodes and edges.
 * <li> <b><code>alphaLengthFile</code></b> - computation of alpha values for semi-Markov states , includes all segments
 * <li> <b><code>expectFile</code></b> - computation of expected values for each Markov feature 
 * <li> <b><code>expectLengthFile</code></b> - computation of expected values for each semi-Markov feature  
 * <li> <b><code>nodeMarginalFile</code></b> - computation of marginal probability of each state at each position 
 * </ul>
 * */
public class CleanMaximumLikelihoodSemiMarkovGradient implements CRFObjectiveFunctionGradient {
	static final Log log = LogFactory.getLog(CleanMaximumLikelihoodSemiMarkovGradient.class);
	public static final boolean debug = log.isDebugEnabled();
	public static final double ASSERTION_TOLERANCE = 0.0001;
	
	public static final int NORM_FACTOR = 50;
	public static final double NORM_MIN = Math.exp(-NORM_FACTOR);
	public static final double NORM_MAX = Math.exp(NORM_FACTOR);

	final LogFiles logs = new LogFiles();
	
	SolverSetup modelInfo;
	CacheProcessor cacheProcessor;
	FeatureEvaluation[] evals;
	LengthFeatureEvaluation[][] lengthEvals;
	boolean[] invalidTransitions;
	
	// / Cache feature information
	// / Cached value of the Mi matrix for all of the features present at every position
	// / Mi is stored as a sparse matrix
	short maxLookback;
	StatePotentials[] statesWithLookback;
	StatePotentials[] statesWithoutLookback;
	int iter = 0;
	double[][] alphas;
	int[] alphaNorms;
	double[] starterAlpha;
	int nSemiMarkovStates;

	// At any given point, lookbackBuffer.get(x) returns the information about a lookback of x. Lookbacks start at 0.
	RecyclingBuffer<LookbackBuffer> lookbackBuffer;
	LookbackBuffer nextBuffer;

	double[] lambda;
	double logZ;
	int zNorm;
	double zInv;
	double[] expects;

	AlphaLengthFeatureProcessor alphaProcessor; 
	BetaLengthFeatureProcessor betaProcessor;
	
	// We publish feature sums 
	private double[] featureSums;
		
	public void setTrainingData(ModelManager fm, List<? extends TrainingSequence<?>> data) {
		cacheProcessor.setTrainingData(fm, data);

		modelInfo = cacheProcessor.getSolverSetup();
		Assert.a(modelInfo.maxStateLengths != null, "Maximum state lengths not set.");
		Assert.a(modelInfo.maxStateLengths.length == modelInfo.nStates, "Maximum state lengths array was length ("+modelInfo.maxStateLengths.length+").  Must have one entry for each state "+modelInfo.nStates+")");
		evals = cacheProcessor.getFeatureEvaluations();
		lengthEvals = cacheProcessor.getLengthFeatureEvaluations();
		invalidTransitions = cacheProcessor.getInvalidTransitions();

		// Create local references to cache values
		maxLookback = modelInfo.maxLookback;
		statesWithLookback = modelInfo.statesWithLookback;
		statesWithoutLookback = modelInfo.statesWithoutLookback;
		nSemiMarkovStates = modelInfo.statesWithLookback.length;

		// Initialize betas (for use later, in the gradient computation)
		alphas = new double[modelInfo.longestSeq][modelInfo.nStates];
		alphaNorms = new int[modelInfo.longestSeq];
		expects = new double[modelInfo.nFeatures];

		LookbackBuffer[] bufferContents = new LookbackBuffer[maxLookback+3];
		for(int i = 0; i<maxLookback+3; ++i) {
			bufferContents[i] = new LookbackBuffer(modelInfo.nStates, modelInfo.nTransitions);
		}
		lookbackBuffer = new RecyclingBuffer<LookbackBuffer>(bufferContents);
		nextBuffer = new LookbackBuffer(modelInfo.nStates, modelInfo.nTransitions);
		
		alphaProcessor = new AlphaLengthFeatureProcessor(this); 
		betaProcessor = new BetaLengthFeatureProcessor(this); 

		starterAlpha = new double[modelInfo.nStates];
	}

	public double apply(double[] param, double[] grad) {
		log.debug(String.format("Beginning It: %d Weights: %s", iter, ColtUtil.format(param)));
		logs.open();
		lambda = param;
		Arrays.fill(grad, 0);
		double totalZ = 0.0;
		double result = 0.0;

		try {
			// Iterate through sequences
			Arrays.fill(expects, 0);
			for (int i = 0; i < modelInfo.nSeqs; ++i) {
				int len = modelInfo.seqOffsets[i + 1] - modelInfo.seqOffsets[i];

				alphaAndBetaPass(i, len);
				
				// Update for the next sequence
				totalZ += logZ;
			}
			
			// sum_j lambda_j F_j(xk, yk)
			double[] featureSums = cacheProcessor.getFeatureSums();
			this.featureSums = featureSums;
			for (int j = 0; j < modelInfo.nFeatures; ++j) {
				result += featureSums[j] * param[j];
				grad[j] = featureSums[j] - expects[j];
			}
			log.debug("Path Value: "+result+" Norm: "+totalZ);
			result -= totalZ;
			if (log.isInfoEnabled()) {
				log.info(String.format("It: %d L=%e, LL=%f, norm(grad): %f Sums: %s Expects: %s Weights: %s Grad: %s", iter, exp(result), result,
						ColtUtil.norm(grad), ColtUtil.format(featureSums), ColtUtil.format(expects), ColtUtil.format(param), ColtUtil.format(grad)));
			}
			Assert.a(exp(result) <= 1.0, "Likelihood is greater than 1.");

			// Normalize by the length of the sequence
			result = result/modelInfo.totalPositions; 
			for(int i=0; i<grad.length; ++i) {
				grad[i] = grad[i]/modelInfo.totalPositions;
			}

			iter += 1;
		}
		finally {
			logs.close();
		}
		return result;
	}

	public void clean() {
	}
		
	void alphaAndBetaPass(int i, int len) {
		// Work forwards, computing alphas
		alphaProcessor.computeAlpha(i, len);
	
		// Since the final beta array is all ones, we can sum the alphas to get the Z
		double sum = 0.0;
		for (double val : alphas[len - 1]) {
			sum += val;
		}
	
		logZ = log(sum) + NORM_FACTOR * (alphaNorms[len - 1]);
		zNorm = ((int) logZ) / NORM_FACTOR;
		zInv = exp(zNorm * NORM_FACTOR - logZ);
		log.debug("Seq: "+i+" Z: "+printNorm(1/zInv, zNorm));
	
		// Work backwards, computing betas and expectations.
		betaProcessor.computeBetasAndExpectations(i, len);
	
		//if(debug) {
		//	logFeatureSums(i);
		//}
	}

	void logBuf() {
		int l = lookbackBuffer.length;
		String s = "";
		for(int i=0; i < l; ++i) {
			s += lookbackBuffer.get(i).pos + " ";
		}
		log.info(s);
	}
	
	void logBufBeta() {
		int l = lookbackBuffer.length;
		String s = "";
		for(int i=0; i < l; ++i) {
			s += ColtUtil.format(lookbackBuffer.get(i).beta) + " ";
		}
		log.info(s);
	}
	
	/**
	 * Computes an unexponentiated mi matrix and updates stable states. Used to create caches for lookback searches.
	 * 
	 */
	void cacheMi(int seqNum, double[] mi, double[] prevStable, double [] newStable, int miPos) {
		if (miPos < 0)
			return;
		//calcMi(mi, overallPosition, cacheStart, cacheStop, false);
		calcMi(mi, seqNum, miPos, false);
		// Go through the mi matrix and for all states with length dependence compute stable values for self transitions
		for (int i = 0; i < modelInfo.nStates; ++i) {
			if (modelInfo.maxStateLengths[i] > 1) {
				// These are all log values so we add them
				newStable[i] = prevStable[i];
				double trans = mi[modelInfo.selfTransitions[i]];
				if(!Double.isInfinite(trans)) {
					//log.debug("Pos: "+miPos+" State: "+i+" Trans: "+trans+" Total: "+(newStable[i]+trans));
					newStable[i] += trans;
				}
			}
		}
	}
	
	/**
	 * This is one of the most time critical parts of the entire solver. The goal is to update the transition matrix.
	 * This function makes a lot of assumptions in order to maximize performance.
	 * 
	 * To maximize performance, we want to make one pass through the Mi matrix, setting each entry to its correct value.
	 * The value for each entry is the exponent of the sum of the weighted feature values of the edge for that entry and
	 * its corresponding node. The entry s0,s1 consists of the s0,s1 edge and the s1 node.
	 */
	void calcMi(double[] mi, int seq, int pos, boolean doExp) {
		cacheProcessor.evaluatePosition(seq, pos);
		double nodeVal = Double.NaN;
		int overallPosition = modelInfo.seqOffsets[seq]+pos;
		int invalidIndex = overallPosition*modelInfo.nPotentials;
		for(short potential : modelInfo.orderedPotentials) {
			boolean invalid = invalidTransitions[invalidIndex + potential];
			double features = invalid ? Double.NEGATIVE_INFINITY : 0.0;

			// Add up all features for this potential.
			FeatureEvaluation potEvals = evals[potential];
			short[] indices = potEvals.index;
			float[] vals = potEvals.value;
			int i = 0;
			short index = indices[i];
			while(index >= 0) {
				// An invalid potential is indicated by a feature value of Short.MAX_VALUE
				features += vals[i]*lambda[index]; 
				index = indices[++i];
			}
			if(index == Short.MIN_VALUE) {
				features = Double.NEGATIVE_INFINITY; 
			}
				
			if(potential < modelInfo.nStates) {
				nodeVal = features;
			}
			else {
				//log.debug(String.format("Mi[%d, %d] = %f, adding in %f to get %f", feat.yprev(), feat.y(), val, feat.value(), val*exp(feat.value()*param[feat.index()])));
				int transition = potential - modelInfo.nStates;
				double val = features + nodeVal;
				if(doExp)
					val = exp(val);
				mi[transition] = val;
			}
		}		
	}

	/**
	 * Given a vector with an existing normalization factor, convert it to a new normalization factor by scaling the
	 * entries.
	 */
	static final void renormalize(final double[] vec, final int currentNorm, final int newNorm) {
		// Instead of dividing by the different (new-current), we reverse the subtraction to negate the exponent and
		// then multiply.
		double factor = exp(NORM_FACTOR * (currentNorm - newNorm));
		//log.info(factor);
		//log.info(ColtUtil.format(vec));
		int len = vec.length;
		for (int i = 0; i < len; ++i) {
			if(vec[i] != 0.0)
				vec[i] *= factor;
		}
		//log.info(ColtUtil.format(vec));
	}

	/** Given a vector, computes a normalization factor for the entries and scales them according to that factor. */
	static final int normalize(final double[] vec) {
		double sum = 0.0;
		for(double val : vec) {
			sum += val;
		}
		if(sum == 0.0 || (sum > NORM_MIN && sum < NORM_MAX)) {
			// No normalization required, our vector is in range.
			return 0;
		}
		if(debug)
			Assert.a(!Double.isNaN(sum));

		//log.info("performing normalization");
		double val = log(sum);
		int norm = (int) val / NORM_FACTOR;
		val = exp(NORM_FACTOR * norm);
		int len = vec.length;
		for (int i = 0; i < len; ++i) {
			vec[i] /= val;
		}
		return norm;
	}

	static final double exp(final double val) {
		return Math.exp(val);
	}

	static final double log(final double val) {
		return Math.log(val);
	}

	final void logFeatureSums(int seqNum) {
		double[][] seqFeatureSums = cacheProcessor.getSequenceFeatureSums();
		if(seqFeatureSums != null) {
			double seqResult = 0.0;
			for (int j = 0; j < modelInfo.nFeatures; ++j) {
				seqResult += seqFeatureSums[seqNum][j] * lambda[j];
			}
			log.debug(String.format("Seq: %d L: %g LL: %f Training path: %f Z: %f", seqNum, exp(seqResult-logZ), seqResult-logZ, seqResult, logZ));
			Assert.a(exp(seqResult-logZ) < 1.0);
		}
	}

	public static final String printNorm(final double value, final int norm) {
		if( value == 0.0)
			return "0 ("+norm+")";
		if( Double.isNaN(value))
			return "NaN ("+norm+")";
		int exponent = (int) log(value);

		double eValue = value/exp(exponent);
		if(Double.isNaN(eValue)) {
			return String.format("NaN(%e n:%d)", value, norm);
		}
		//return String.format("%e(%d) %fe%d", value, norm, eValue, exponent+norm*NORM_FACTOR);
		return String.format("%fe%d", eValue, exponent+norm*NORM_FACTOR);
	}

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

	public String getAlphaLengthFile() {
		return logs.alphaLengthFile;
	}

	public void setAlphaLengthFile(String alphaLengthFile) {
		logs.alphaLengthFile = alphaLengthFile;
	}

	public String getAlphaFile() {
		return logs.alphaFile;
	}

	public void setAlphaFile(String alphaFile) {
		logs.alphaFile = alphaFile;
	}

	public String getExpectFile() {
		return logs.expectFile;
	}

	public void setExpectFile(String expectFile) {
		logs.expectFile = expectFile;
	}

	public String getExpectLengthFile() {
		return logs.expectLengthFile;
	}

	public void setExpectLengthFile(String expectLengthFile) {
		logs.expectLengthFile = expectLengthFile;
	}

	public String getNodeMarginalFile() {
		return logs.nodeMarginalFile;
	}

	public void setNodeMarginalFile(String nodeMarginalFile) {
		logs.nodeMarginalFile = nodeMarginalFile;
	}

	public String getBetaLengthFile() {
		return logs.betaLengthFile;
	}

	public void setBetaLengthFile(String betaLengthFile) {
		logs.betaLengthFile = betaLengthFile;
	}
}
