package calhoun.analysis.crf.solver.check;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.CRFObjectiveFunctionGradient;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.analysis.crf.solver.CacheProcessor;
import calhoun.analysis.crf.solver.RecyclingBuffer;
import calhoun.util.Assert;
import calhoun.util.ColtUtil;
import calhoun.util.DenseIntMatrix2D;

/**
 * The general normalization scheme works as follows. When updating alpha values in the forward pass we compute segments
 * of length 1 first and then work backwards.
 * 
 * Instead of always normalizing to 1 we discretize the normalization. We choose an arbitrary normalization factor w,
 * such as 50. The normalization factor at any position is then an integer v, and all entries at that position are
 * alpha[y]*e^(v*w).
 * 
 * The normalization can be computed at any position from 1) Elements of the alpha array are summed s 2) v = log(s)/w.
 * By integer division v will always be an appropriate normalizer. It may be positive or negative. 3) All elements of
 * the array are divided by e^(v*w)
 * 
 * To add values
 */
public class CachedSemiCRFGradient implements CRFObjectiveFunctionGradient {
	private static final Log log = LogFactory.getLog(CachedSemiCRFGradient.class);
	boolean debug = log.isDebugEnabled();
	static final int NORM_FACTOR = 10;
	static final double NORM_MIN = Math.exp(-NORM_FACTOR);
	static final double NORM_MAX = Math.exp(NORM_FACTOR);
	// / Cache feature information
	short[] id;
	byte[] potentialIx;
	float[] val;
	// / Cached value of the Mi matrix for all of the features present at every position
	// / Mi is stored as a sparse matrix
	int nTransitions;
	short[] transitionFrom;
	short[] transitionTo;
	short[] orderedPotentials;
	boolean[] invalidTransitions;
	DenseIntMatrix2D transitionIndex;
	// / Cached values of the sums of each feature value through the whole training set.
	double[] featureSums;
	// / Index into the feature arrays of the first feature for each postion of each sequence.
	int[] starts;
	// / Index into the starts array of the first position of each sequence.
	int[] seqOffsets;
	// / Number of sequences in the training data set
	int nSeqs;
	// / Number of constant features.
	int nConstantFeatures;
	short maxStateLengths[];
	short maxLookback;
	int[] lengthStarts;
	short[] lookbacks;
	byte[] lengthPotentials;
	short[] lengthIndexes;
	float[] lengthVals;
	CacheProcessor.StatePotentials[] statesWithLookback;
	CacheProcessor.StatePotentials[] statesWithoutLookback;
	List<? extends TrainingSequence<?>> data;
	ModelManager fm;
	int nFeatures;
	int nStates;
	int nPotentials;
	int iter = 0;
	double[][] alphas;
	int[] alphaNorms;
	double[] starterAlpha;

	/**
	 * This object holds information about previous positions during the computation of betas and expectations. This
	 * allows us to quickly access data about previous positions. These objects are kept in a recycling buffer that
	 * keeps one buffer for each possible lookback.
	 * 
	 * One tricky aspect of this is that the details change slightly between the forward and backwards pass.  On the forward
	 * pass, the lookback contains the information in the normal way.  In the backwards pass, stable states and transitions are 
	 * shifted back one base compared to the betas.
	 */
	class LookbackBuffer {
		int pos;
		
		// The mi matrix for transitioning from pos-lookback-1 to pos-lookback
		double[] mi = new double[nPotentials];

		// The weighted sum of feature values for staying in this position from the beginning to pos-lookback
		double[] stableState = new double[nStates];

		// Initial values of the beta vector for somelength dependent states.
		double[] beta = new double[nStates];
		
		// Norm of prevBeta.
		int betaNorm;
		
		// Stores the probability of all segments begining at this position using this transition.
		double[] transitionProb = new double[nTransitions];

		/** mi and stableStates are cleared as new values are entered. This fixes the others */
		void clear()
		{
			pos = -1;
			Arrays.fill(beta, 0.0);
			betaNorm = 0;
			Arrays.fill(transitionProb, 0.0);
		}
	}

	// At any given point, lookbackBuffer.get(x) returns the information about a lookback of x. Lookbacks start at 0.
	RecyclingBuffer<LookbackBuffer> lookbackBuffer;
	LookbackBuffer nextBuffer;

	double[] lambda;
	double[] constMi;
	double logZ;
	int zNorm;
	double zInv;
	double[] expects;

	AlphaLengthFeatureProcessor alphaProcessor; 
	BetaLengthFeatureProcessor betaProcessor; 
	
	double exp(double argVal) {
		return Math.exp(argVal);
	}

	double log(double argVal) {
		return Math.log(argVal);
	}

	boolean allPaths;
	public CachedSemiCRFGradient(short[] maxStateLengths, boolean allPaths) {
		this.maxStateLengths = maxStateLengths;
		this.allPaths = allPaths;
	}

	public void clean() {
	}
	
	public void setTrainingData(ModelManager fm, List<? extends TrainingSequence<?>> data) {
		this.fm = fm;
		this.data = data;
		nFeatures = fm.getNumFeatures();
		nStates = fm.getNumStates();
		nSeqs = data.size();
		
		FeatureCacheLength cache = new FeatureCacheLength(fm, data, allPaths, maxStateLengths, new short[maxStateLengths.length], true);
		// Create local references to cache values
		orderedPotentials = cache.orderedPotentials;
		id = cache.id;
		potentialIx = cache.potentialIx;
		val = cache.val;
		transitionFrom = cache.transitionFrom;
		transitionTo = cache.transitionTo;
		nTransitions = cache.nTransitions;
		transitionIndex = cache.transitionIndex;
		featureSums = cache.featureSums;
		starts = cache.starts;
		seqOffsets = cache.seqOffsets;
		invalidTransitions = cache.invalidTransitions;
		nConstantFeatures = cache.numConstantFeatures;
		nPotentials = cache.nPotentials;
		maxLookback = cache.maxLookback;
		lengthStarts = cache.lengthStarts;
		lookbacks = cache.lookbacks;
		lengthPotentials = cache.lengthPotentials;
		lengthIndexes = cache.lengthIndexes;
		lengthVals = cache.lengthVals;
		statesWithLookback = cache.statesWithLookback;
		statesWithoutLookback = cache.statesWithoutLookback;

		// Initialize betas (for use later, in the gradient computation)
		alphas = new double[cache.longestSeq][nStates];
		alphaNorms = new int[cache.longestSeq];
		constMi = new double[nTransitions];
		expects = new double[nFeatures];

		LookbackBuffer[] bufferContents = new LookbackBuffer[maxLookback+3];
		for(int i = 0; i<maxLookback+3; ++i) {
			bufferContents[i] = new LookbackBuffer();
		}
		lookbackBuffer = new RecyclingBuffer<LookbackBuffer>(bufferContents);
		nextBuffer = new LookbackBuffer();
		
		alphaProcessor = new AlphaLengthFeatureProcessor(); 
		betaProcessor = new BetaLengthFeatureProcessor(); 

		starterAlpha = new double[nStates];
	}

	public double apply(double[] param, double[] grad) {
		lambda = param;
		Arrays.fill(grad, 0);
		double result = 0.0;

		// Calculate the constant Mi matrix
		// CalcMi produces an mi matrix, so we copy them to the constMi matrix.
		Arrays.fill(constMi, 0.0);
		calcMi(constMi, -1, 0, starts[0], false);

		// Iterate through sequences
		Arrays.fill(expects, 0);
		int seqStart = 0;
		for (int i = 0; i < nSeqs; ++i) {
			int len = seqOffsets[i + 1] - seqOffsets[i];
			// Work forwards, computing alphas
			alphaProcessor.computeAlpha(seqStart, len);

			// Since the final beta array is all ones, we can sum the alphas to get the Z
			double sum = 0.0;
			for (double localVal : alphas[len - 1]) {
				sum += localVal;
			}

			logZ = log(sum) + NORM_FACTOR * (alphaNorms[len - 1]);
			zNorm = ((int) logZ) / NORM_FACTOR;
			zInv = exp(zNorm * NORM_FACTOR - logZ);
			//log.info(String.format("Z: %f (Norm: %d)", 1/zInv, zNorm));

			// Work backwards, computing betas and expectations.
			betaProcessor.computeBetasAndExpectations(seqStart, len);

			// Update for the next sequence
			result -= logZ;
			seqStart += len;
		}
		// sum_j lambda_j F_j(xk, yk)
		for (int j = 0; j < nFeatures; ++j) {
			result += featureSums[j] * param[j];
			grad[j] = featureSums[j] - expects[j];
		}
		if (log.isInfoEnabled()) {
			log.info(String.format("It: %d L=%e, LL=%f, norm(grad): %f Sums: %s Expects: %s Weights: %s Grad: %s", iter, exp(result), result,
					ColtUtil.norm(grad), ColtUtil.format(featureSums), ColtUtil.format(expects), ColtUtil.format(param), ColtUtil.format(grad)));
			// log.info(String.format("It: %d L=%e, LL=%f, norm(grad): %f Weights: %s Grad: %s", iter, exp(result),
			// result, ColtUtil.norm(grad), ColtUtil.format(param), ColtUtil.format(grad)));
		}
		iter += 1;
		int totalPositions = seqOffsets[seqOffsets.length-1];
		result = result/totalPositions; 
		for(int i=0; i<grad.length; ++i) {
			grad[i] = grad[i]/totalPositions;
		}
		return result;
	}

	class AlphaLengthFeatureProcessor extends LengthFeatureProcessor {
		int pos;
		double[] alpha;
		int alphaNorm;
		double[] stableState;
		
		/**
		 * In the forward pass we compute alpha values and expections. This is simpler than the backwards pass because
		 * the cache is set up for us so that we can always look at one position at a time. We have to cache previous
		 * values but we never have to look ahead.
		 */
		void computeAlpha(int seqStart, int len) {
			Arrays.fill(alphaNorms, 0);
			Arrays.fill(starterAlpha, 0.0);
			int cacheStart = starts[seqStart];
			double[] prevAlpha = null;

			for(pos = 0; pos < len; ++pos) {
				int overallPosition = seqStart + pos;
				int cacheStop = starts[overallPosition+1];
	
				prevAlpha = alpha;
				alpha = alphas[pos];
				Arrays.fill(alpha, 0.0);
				if (pos == 0) {
					alphaNorm = 0;
					calcStartAlpha(alpha, overallPosition, cacheStart, cacheStop);
					// Put an empty entry in the lookback so the first base has 0's initialized.
					Arrays.fill(nextBuffer.stableState, 0.0);
				} else {
					cacheMi(nextBuffer.mi, stableState, nextBuffer.stableState, seqStart, pos, cacheStart, cacheStop);
					regularAlphaUpdate(pos, nextBuffer.mi, prevAlpha, alpha);
				}

				// Add the lookback into the array
				stableState = nextBuffer.stableState;
				nextBuffer = lookbackBuffer.addFirst(nextBuffer);
					
				// Now we need to loop through for the length dependent cache
				alphaProcessor.lengthCache(overallPosition);
	
				alphaNorm += normalize(alpha);
				alphaNorms[pos] = alphaNorm;
				//log.info(String.format("alphaNorm[%d]=%d", pos, alphaNorm));
				cacheStart = cacheStop;
			}
		}

		/**
		 * Updates the alpha vector for non-length dependent states. We don't have to worry about normalization here
		 * because the regular alpha update is done before the length dependent, so these are the first values that will
		 * be set.
		 */
		private void regularAlphaUpdate(int argPos, double[] mi, double[] lastAlpha, double[] newAlpha) {
			double nodeVal = 0.0;
			int lastState = -1;
			boolean lengthNode = false;
			for (short potential : orderedPotentials) {
				if (potential < nStates) {
					if (lastState != -1) {
						newAlpha[lastState] = nodeVal;
					}
					lastState = potential;
					nodeVal = 0.0;
					lengthNode = maxStateLengths[potential] > 1;
				} else {
					if (!lengthNode) {
						int trans = potential - nStates;
						int from = transitionFrom[trans];
						//log.debug(String.format("alpha[%d][%d] = %f = %f + alpha[%d][%d] %f * %f (%f)", pos, lastState, nodeVal + lastAlpha[from] * exp(mi[trans]), nodeVal, pos-1, from, lastAlpha[from], exp(mi[trans]), mi[trans]));
						nodeVal += lastAlpha[from] * exp(mi[trans]);
					}
				}
			}
			newAlpha[lastState] = nodeVal;
		}

		/** Updates an alpha entry with a weighted sum of features values for a given potential */
		@Override void potential(short lookback, byte toNode, int trans, float potentialValue, int nodeStart, int nodeStop, int edgeStart, int edgeStop) {
			//log.debug(String.format("Alpha Length cache: Pos: %d Lb: %d State: %d Trans: %d Pot: %f Node: %d-%d Edge: %d-%d", pos, lookback, toNode, trans, potentialValue, nodeStart, nodeStop, edgeStart, edgeStop));
			/*
			 * Updates an existing alpha by adding in: potentialValue - The value of any length-dependent features for
			 * this node f(y, i, d) and edge f(y', y, i, d) stableValue - The value of the non-length dependent node
			 * features summed across the length of this segment mis - The value of the non-length dependent transition
			 * from the previous node to this one f(y', y, i-d)
			 */
			int prevPos = pos - lookback - 1;
			LookbackBuffer buffer = lookbackBuffer.get(lookback);

			double stableValue = stableState[toNode] - buffer.stableState[toNode];
			double expVal = potentialValue + stableValue;
			double prevAlpha = 1.0;
			if(prevPos >= 0) {
				int fromNode = transitionFrom[trans];
				expVal += buffer.mi[trans] + alphaNorms[prevPos] * NORM_FACTOR;
				prevAlpha = alphas[prevPos][fromNode];
			}
			else {
				expVal += starterAlpha[toNode];
			}
			int norm = ((int) expVal) / NORM_FACTOR;
			expVal -= norm * NORM_FACTOR;

			if (Math.abs(norm) > Math.abs(alphaNorm)) {
				renormalize(alpha, alphaNorm, norm);
				//log.debug("Renormalized: "+ColtUtil.format(alpha));
				alphaNorm = norm;
			} else if (Math.abs(norm) < Math.abs(alphaNorm)) {
				expVal += NORM_FACTOR * (norm - alphaNorm);
			}
			double update = exp(expVal) * prevAlpha;
			alpha[toNode] += update;
			//log.debug(String.format("alpha[%d][%d] = %f (%d) = %f + %f (alpha[%d][%d]) * %f (Pot: %f Stab: %f Trans: %f Norm: %f)", pos, toNode, alpha[toNode], alphaNorm, alpha[toNode]-update, prevAlpha, pos-lookback-1, trans == -1 ? -1 : transitionFrom[trans], exp(expVal), potentialValue, stableValue, trans == -1 ? starterAlpha[toNode] : buffer.mi[trans], trans == -1 ? 0.0 : alphaNorms[prevPos] * NORM_FACTOR));
		}

		/**
		 * A specialized version of calcMi for the first position in a sequence. Has the special property that constant
		 * edge features are not included.
		 */
		void calcStartAlpha(double[] alpha1, int overallPosition, int posCurrent, int posStop) {
			// Constant features
			int constCurrent = 0;
			byte constPotential = -1;
			double constVal = Double.NaN;
			// Positional features
			byte posPotential = -1;
			double posVal = Double.NaN;
			int invalidIndex = overallPosition * nPotentials;
			for (short potential : orderedPotentials) {
				boolean invalid = invalidTransitions[invalidIndex + potential];
				double features = invalid ? Double.NEGATIVE_INFINITY : 0.0;
				// Find constant features for this potential. Ignore edge features.
				while (constPotential == -1 || constPotential == potential) {
					if (constPotential != -1 && potential < nStates) {
						features += constVal;
						Assert.a(!Double.isNaN(features));
					}
					if (constCurrent < nConstantFeatures) {
						constVal = val[constCurrent] * lambda[id[constCurrent]];
						constPotential = potentialIx[constCurrent];
						++constCurrent;
					} else {
						break;
					}
				}
				if (potential < nStates) {
					// Add in cached features, but only for non-length dependent nodes
					while (posPotential == -1 || posPotential == potential) {
						if (posPotential != -1)
							features += posVal;
						if (posCurrent < posStop) {
							posVal = val[posCurrent] * lambda[id[posCurrent]];
							posPotential = potentialIx[posCurrent];
							++posCurrent;
						} else {
							break;
						}
					}
					if(maxStateLengths[potential]> 1) {
						starterAlpha[potential] = features;
					}
					else {
						alpha1[potential] = exp(features);
					}
				}
			}
			// verify that we actually checked all the cache entries.
			Assert.a(constCurrent == nConstantFeatures);
			Assert.a(posCurrent == posStop);
		}
	}

	/** Updates a beta entry with a potential from a length dependent feature */
	class BetaLengthFeatureProcessor extends LengthFeatureProcessor {
		int lengthPos;
		double[] lengthStable;
		int miPos;
		double[] stableState;

		double[] beta;
		int betaNorm;

		LookbackBuffer posLookback; 
		LookbackBuffer prevLookback; 
		
		// This vector tracks the marginal probability of each node, accounting for length dependence
		double[] nodeProb = new double[nStates];
		double[] newNodeProb = new double[nStates];

		// This vector tracks the marginal probability of each edge, accounting for length dependence
		double[] edgeProb = new double[nTransitions];
		
		/**
		 * Computing the betas is a bit tricky. The problem is that our length based cache associates each explicit
		 * length feature with the position at which it ends. For the betas, we are going backwards, and so we need to
		 * look at the positions where the length based feature starts.
		 * 
		 * The way we handle this is by evaluating the length based features as they arise and incrementally adding
		 * their potential values into the beta positions as we go along. Therefore we aren't always filling one beta
		 * matrix at a time, we fill any beta entry that starts a state which ends the current position as we move
		 * backwards through the sequence.
		 * 
		 * In order to do this we need to keep mi matrices and arrays of the potentials of staying in a stable state for
		 * the entire lookback period. These matrices let us compute betas anywhere we want.
		 * 
		 * Unfortunately this process needs more exponentiation than the existing process.
		 * 
		 * Normalization works as follows: First, we need to find the right normalization factor for the value we are
		 * adding. 1) We add all of the exponents we will need. 2) We integer divide by NORM_FACTOR to get the
		 * normalization constant for our beta update 3) The constant*NORM_FACTOR is subtracted from the total exponent
		 * to get the unnormalized exponent 4) The unnormalized exponent is exponentiated and multiplied by the
		 * following beta value 5) We compare the norm constant to the current constant for the beta we are about to
		 * update 5.1) If they are equal we are fine, the values can be added 5.2) If Our new norm constant has larger
		 * abs(), we renormalize the beta vector by dividing all elements by e^(NORM_FACTOR*diff) 5.3) If our new norm
		 * constant has smaller abs(), we renormalize it but dividing by e^(FACTOR*diff) 6) When we finally update the
		 * beta vector at the end we redo the normalization anyway.
		 * 
		 * To optimize this for fastest execution it ends up being one big function. Putting function calls in the inner
		 * loops really slowed things down. I think the optimizer can't do as good a job if there are function calls
		 * there.
		 */
		void computeBetasAndExpectations(int seqStart, int len) {
			// Start by initializing an array of unexponentiated mi and stable state matrices for all of the possible
			// lookback positions from the last base
			// Also initialize the array of previous beta vectors
			int lastInitPos = len - 2 - maxLookback;
			int cacheStop = starts[seqStart + len];

			// posLookback holds the lookback buffer for the position currently specified by pos
			posLookback = nextBuffer;
			// prevLookback holds the lookback buffer for the position pos+1
			prevLookback = null;

			Arrays.fill(nodeProb, 0.0);
			Arrays.fill(newNodeProb, 0.0);
			
			// Now work through the sequence backwards
			miPos = len-1;
			nextBuffer.clear();
			for (int pos = len - 1; pos >= 0; --pos) {
				// First, update the lookback, which caches mi and stable values if necessary
				while(miPos >= 0 && miPos >= lastInitPos) {
					if(miPos == len-1) {
						// Initialize stable states
						// Need to prime with a stable state of 0 since we don't have the init that we do for alphas.
						Arrays.fill(nextBuffer.stableState, 0);
					}
					else {
						// Update stable states given the previous transition
						int cacheStart = starts[seqStart + miPos+1];
						nextBuffer.clear();
						cacheMi(nextBuffer.mi, stableState, nextBuffer.stableState, seqStart, miPos+1, cacheStart, cacheStop);
						cacheStop = cacheStart;
					}
					nextBuffer.pos = miPos;
					stableState = nextBuffer.stableState;
					nextBuffer = lookbackBuffer.addFirst(nextBuffer);
					nextBuffer.clear();
					--miPos;
				}
				//logBuf();
				//log.debug("Pos: "+pos+" MiPos: "+miPos);
				posLookback = lookbackBuffer.get(pos - miPos -1);

				// Check that posLookback is correct
				Assert.a(posLookback.pos == pos, "Wrong lookback buffer: was ", posLookback.pos, " should be ", pos);

				// Update betas (if necessary)
				if(prevLookback == null) {
					// For the first iteration (last position in the sequence)
					Assert.a(pos == len-1);

					// Initialize first betas
					Arrays.fill(posLookback.beta, 1);
					posLookback.betaNorm = 0;

					// Initialize node marginals
					double nodeNorm = exp((alphaNorms[pos] - zNorm) * NORM_FACTOR) * zInv;
					for(int i=0; i<nStates; ++i) {
						nodeProb[i] = nodeNorm * alphas[pos][i];
					}
				}
				else {
					Assert.a(prevLookback.pos == pos+1);
					posLookback.betaNorm = regularBetaUpdate(pos+1, posLookback.beta, posLookback.betaNorm, prevLookback.beta, prevLookback.betaNorm, prevLookback.transitionProb, posLookback.mi);
					//log.info(String.format("Pos: %d, Mipos: %d LbPos: %d", pos, miPos, posLookback.pos));
					
					// Now we need to loop through for the length dependent cache
					beta = prevLookback.beta;
					betaNorm = prevLookback.betaNorm;
					lengthStable = prevLookback.stableState;
					lengthPos = pos + 1;
					betaProcessor.lengthCache(seqStart + lengthPos);

					System.arraycopy(nodeProb, 0, newNodeProb, 0, nStates);
					
					// Now calculate edge marginals for staying in the same state
					// Take the node probabilities, and subtract off each transition
					for(CacheProcessor.StatePotentials lb: statesWithLookback) {
						int state = lb.state;
						int index = transitionIndex.getQuick(state, state);
						double transProb = nodeProb[state];
						for(byte pot : lb.potentials) {
							double lbTrans = posLookback.transitionProb[pot - nStates];
							transProb -= lbTrans;
							newNodeProb[state] -= lbTrans;
						}
						Assert.a(posLookback.transitionProb[index] == 0.0);
						posLookback.transitionProb[index] = transProb;
					}
				}
				
				// As a check, we verify that the node marginals sum to one for each position.
				double sum =0.0;
				for(double x : nodeProb) {
					sum += x;
				}
				if(Math.abs(1.0-sum) > 0.001) {
					Assert.a(false, "Node marginals don't sum to 1: ", ColtUtil.format(nodeProb));
				}

				// Verify that edge marginals sum to 1.
				if(prevLookback != null) {
					sum =0.0;
					for(double x : posLookback.transitionProb) {
						sum += x;
					}
					for(double x : edgeProb) {
						sum += x;
					}
					if(Math.abs(1.0-sum) > 0.001) {
						Assert.a(false, "Edge marginals don't sum to 1: ", ColtUtil.format(edgeProb), ColtUtil.format(posLookback.transitionProb));
					}

					// At this point the previous beta values are all updated.
					// Update expectations for the transitions we just calculated.
					updateExpectations(seqStart, pos+1, posLookback.transitionProb);
				}

				// Update the node probabilities to remove the length transitions
				double[] temp = nodeProb;
				nodeProb = newNodeProb;
				newNodeProb = temp;
				
				if (debug) {
					if ((seqStart == 0) && (pos < 2 || pos >= len - 2)) {
						log.debug(String.format("Pos: %d expects: %s alphas: %s (norm %d) betas: %s (norm %d) MiPos: %d", pos, ColtUtil.format(expects), ColtUtil
								.format(alphas[pos]), alphaNorms[pos], ColtUtil.format(posLookback.beta), posLookback.betaNorm, miPos+1));
					}
				}

				prevLookback = posLookback;
				--lastInitPos;
			}
			// Now update for the first position
			//log.debug(String.format(ColtUtil.format(posLookback.beta)));
			posLookback.betaNorm = regularBetaUpdate(0, null, posLookback.betaNorm, posLookback.beta, posLookback.betaNorm, null, null);
			beta = posLookback.beta;
			betaNorm = posLookback.betaNorm;
			lengthStable = posLookback.stableState;
			lengthPos = 0;
			betaProcessor.lengthCache(seqStart);
			updateExpectations(seqStart, 0, posLookback.transitionProb);
		}

		/**
		 * Does the update of the beta values for all of the regular, non-length dependent states. Also calculates all
		 * of the node and edge probabilities for the non-length nodes and the transitions into them.
		 */
		private int regularBetaUpdate(int pos, double[] newBeta, int newNorm, double[] oldBeta, int oldNorm, double[] transitionProb, double[] mi) {
			// Need to deal with possibly different normalization constants.
			int norm = newNorm;
			double normAdjust = 0.0;
			if (Math.abs(oldNorm) > Math.abs(newNorm)) {
				// We never make the constant smaller, so set the constant for the new factor
				//log.info(String.format("Renormalizing beta[%d] from %d to %d", pos, newNorm, oldNorm));
				renormalize(newBeta, newNorm, oldNorm);
				norm = oldNorm;
				newNorm = oldNorm;
			} else {
				// The case where beta(pos-1) already has a normalization constant larger than beta(pos)
				// We need to adjust all of the updated beta values.
				normAdjust = (oldNorm - newNorm) * NORM_FACTOR;
			}

			//log.info(String.format("Node norm pos: %d e^(alpha: %d + beta: %d - z: %d) * zinv: %f ", pos, alphaNorms[pos], norm, zNorm, zInv));
			double[] nodeAlpha = alphas[pos];
			double nodeNorm = exp((alphaNorms[pos] + oldNorm - zNorm) * NORM_FACTOR) * zInv;
			double[] edgeAlpha = null;
			double edgeNorm = Double.NaN;
			if(pos > 0) {
				edgeAlpha = alphas[pos-1];
				// We add newNorm here because it ends up cancelling with normAdjust when we calc the edge prob.
				edgeNorm = exp((alphaNorms[pos-1] + newNorm - zNorm) * NORM_FACTOR) * zInv;
			}

			for(CacheProcessor.StatePotentials potentials : statesWithoutLookback) {
				byte node = potentials.state;
				double nodePotential = 0.0;
				double betaVal = oldBeta[node];
				//log.info(String.format("NodeMarg[%d][%d] = %f = %f * %f * %f (aN: %d bN: %d zN: %d 1/z: %f)", pos, node, nodeAlpha[node] * betaVal * nodeNorm, nodeAlpha[node], betaVal, nodeNorm, alphaNorms[pos], oldNorm, zNorm, zInv));
				nodeProb[node] = nodeAlpha[node] * betaVal * nodeNorm;
				
				// For regular states, we sum edge probabilities to get node probabilities.
				if(pos > 0) {
					for (short potential : potentials.potentials) {
						int trans = potential - nStates;
						int from = transitionFrom[trans];
						// Mi for beta is not exponentiated, so we do it here.
						double potentialValue = exp(mi[trans] + normAdjust);
						nodePotential += potentialValue;
						//log.info(String.format("Beta[%d][%d] = %f = %f + Beta[%d][%d] %f * %f (Mi: %f)", pos, from, newBeta[from] + potentialValue * betaVal, newBeta[from], pos+1, node, betaVal, potentialValue, mi[trans]));
						newBeta[from] += potentialValue * betaVal;
						edgeProb[trans] = edgeAlpha[from] * potentialValue * betaVal * edgeNorm;
					}
				}
			}

			// Now check to see if this new vector needs normalization again.
			int ret = norm;
			if(newBeta != null) {
				ret += normalize(newBeta);
			}
			return ret;
		}

		/**
		 * Updates an existing beta by adding in: potentialValue - The value of any length-dependent features for this
		 * node f(y, i, d) and edge f(y', y, i, d) stableValue - The value of the non-length dependent node features
		 * summed across the length of this segment mis - The value of the non-length dependent transition from the
		 * previous node to this one f(y', y, i-d)
		 * 
		 * In addition to computing the betas for length dependent features, it updates the probabilities vectors needed
		 * for the feature expectations. These are the marginal at each position or each edge.
		 */
		@Override void potential(short lookback, byte toNode, int trans, float potentialValue, int nodeStart, int nodeStop, int edgeStart, int edgeStop) {
			//log.debug(String.format("Beta Length cache: Pos: %d Lb: %d State: %d Trans: %d Pot: %f Node: %d-%d Edge: %d-%d", lengthPos, lookback, toNode, trans, potentialValue, nodeStart, nodeStop, edgeStart, edgeStop));
			int bufferPos = lengthPos - lookback - 1;
			int lbIndex = bufferPos - miPos - 1;
			//log.info("Mipos "+miPos+" index "+lbIndex);
			LookbackBuffer buffer = trans == -1 ? null : lookbackBuffer.get(lbIndex);
			LookbackBuffer stableBuffer = lookbackBuffer.get(lbIndex+1);
			int prevPos = lengthPos - lookback - 1;
			int fromNode = trans == -1 ? -1 : transitionFrom[trans];

			//Assert.a(lengthPos == 0 || stabBuffer.pos == lengthPos + 1 -lookback, "Bad stab buf ", buffer.pos, " not ", lengthPos+1-lookback);
			double stableValue = stableBuffer.stableState[toNode] - lengthStable[toNode];

			double prevAlpha = 1.0;
			double prevAlphaNorm = 0;
			double transVal;
			if(trans != -1) {
				Assert.a(prevLookback.pos == lengthPos, "Expected ",lengthPos, " was ",prevLookback.pos);
				Assert.a(buffer.pos == (lengthPos-lookback-1), "Expected ",(lengthPos-lookback-1), " was ",buffer.pos);
				transVal = buffer.mi[trans];
				prevAlpha = alphas[prevPos][fromNode];
				prevAlphaNorm = alphaNorms[prevPos];
			}
			else {
				transVal = starterAlpha[toNode];
			}
			double expVal = potentialValue + stableValue + transVal;

			int norm = ((int) expVal) / NORM_FACTOR;
			expVal -= norm * NORM_FACTOR;
			norm += betaNorm;
			
			// In addition to updating the beta array, we need to calculate a probability for this segment so we can
			// correctly calculate feature expectations
			double prob = prevAlpha * beta[toNode]
					* exp(expVal + NORM_FACTOR * (prevAlphaNorm + norm - zNorm)) * zInv;
			Assert.a(!Double.isNaN(prob));
			
			// Now update expectations for all node features for this edge
			int currentFeature = nodeStart;
			while (currentFeature < nodeStop) {
				short index = lengthIndexes[currentFeature];
				if(index != -1) {
					expects[index] += prob * lengthVals[currentFeature];
				}
				++currentFeature;
			}

			currentFeature = edgeStart;
			while (currentFeature < edgeStop) {
				short index = lengthIndexes[currentFeature];
				if(index != -1) {
					expects[index] += prob * lengthVals[currentFeature];
				}
				++currentFeature;
			}

			//log.info(String.format("NodeMarg[%d][%d] = %f = %f + %f (p: %f mi: %f s: %f) * %f alpha[%d][%d] * %f * %f", lengthPos, toNode, nodeProb[toNode]+prob, nodeProb[toNode], exp(expVal), potentialValue, transVal, stableValue, prevAlpha, bufferPos, fromNode, beta[toNode], zInv));
			nodeProb[toNode] += prob;

			// Update the beta values
			// Updates the transition probabilities
			if(trans != -1) {
				if (Math.abs(norm) > Math.abs(buffer.betaNorm)) {
					renormalize(buffer.beta, buffer.betaNorm, norm);
					buffer.betaNorm = norm;
				} else if (Math.abs(norm) < Math.abs(buffer.betaNorm)) {
					expVal -= NORM_FACTOR * (buffer.betaNorm - norm);
				}
				double transPotential = exp(expVal);

				buffer.transitionProb[trans] += prob;
				double update = transPotential * beta[toNode];
				buffer.beta[fromNode] += update;
				//log.info(String.format("Beta[%d][%d] = %f (%d) = %f + %f (beta[%d][%d]) * %f (Pot: %f Stab: %f Trans: %f)", bufferPos, fromNode, buffer.beta[fromNode], buffer.betaNorm, buffer.beta[fromNode]-update, beta[toNode], lengthPos, toNode, transPotential, potentialValue, stableValue, trans == -1 ? 0.0 : buffer.mi[trans]));
			}
		}

		void updateExpectations(int seqStart, int pos, double[] transitionProb) {
			//log.info(String.format("Pos %d: %s", pos, ColtUtil.format(transitionProb)));
			int overallPosition = seqStart + pos;
			int posCurrent = starts[overallPosition];
			int posStop = starts[overallPosition + 1];

			// First compute the expectations for the length-dependent states.
			// Constant features
			int constCurrent = 0;
			short constId = -1;
			byte constPotential = -1;
			double constVal = Double.NaN;

			// Positional features
			short posId = -1;
			byte posPotential = -1;
			double posVal = Double.NaN;

			boolean includeEdges = pos != 0;
			int invalidIndex = overallPosition * nPotentials;

			boolean lengthNode = false;
			for (short potential : orderedPotentials) {
				boolean invalid = invalidTransitions[invalidIndex + potential];

				double prob = Double.NaN;
				if(potential < nStates) {
					lengthNode = maxStateLengths[potential] > 1;
					prob = nodeProb[potential];
				}
				else {
					prob = (lengthNode ? transitionProb : edgeProb)[potential-nStates]; 
					//prob = transitionProb[potential-nStates]; 
				}
				
				// Include constant features for this potential. Need to skip over constant features for invalid
				// potentials
				while (constPotential == -1 || constPotential == potential) {
					if (constPotential != -1 && !invalid && (includeEdges || potential < nStates)) {
						//log.info(String.format("Expect #%d: %f * %f ", constId, prob, constVal));
						expects[constId] += prob * constVal;
					}
					if (constCurrent < nConstantFeatures) {
						constId = id[constCurrent];
						constVal = val[constCurrent];
						constPotential = potentialIx[constCurrent];
						++constCurrent;
					} else {
						break;
					}
				}

				// Include cached features for this potential. The cache should never have features for invalid
				// potentials.
				if (!invalid) {
					while (posId == -1 || posPotential == potential) {
						if (posId != -1) {
							//log.info(String.format("Expect #%d: %f * %f ", posId, prob, posVal));
							expects[posId] += prob * posVal;
						}
						if (posCurrent < posStop) {
							posId = id[posCurrent];
							posVal = val[posCurrent];
							posPotential = potentialIx[posCurrent];
							++posCurrent;
						} else {
							break;
						}
					}
				}
			}
			// verify that we actually checked all the cache entries.
			Assert.a(constCurrent == nConstantFeatures);
			Assert.a(posCurrent == posStop);
		}
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
	void cacheMi(double[] mi, double[] prevStable, double [] newStable, int seqStart, int miPos, int cacheStart, int cacheStop) {
		if (miPos < 0)
			return;
		int overallPosition = seqStart + miPos;
		calcMi(mi, overallPosition, cacheStart, cacheStop, false);
		// Go through the mi matrix and for all states with length dependence compute stable values for self transitions
		for (int i = 0; i < nStates; ++i) {
			if (maxStateLengths[i] > 1) {
				// These are all log values so we add them
				newStable[i] = mi[transitionIndex.getQuick(i, i)] + prevStable[i];
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
	 * 
	 * Because node features are applied to more than 1 entry in the matrix, we use a sorting of all of the features
	 * where each node preceeds all its corresponding edges. This allows us to keep track of only 1 node value at a time
	 * and easily apply it to all its edge features.
	 * 
	 * As we evaluate each potential, we check the cache to see if it is valid at this position and to get any features
	 * values. Note that this function very much depends on the fact that the entries in the cache will be in the
	 * correct order.
	 * 
	 * The other wrinkle is that for features that always occur (constant features), we pull them from the constant mi
	 * array, not from the cache.
	 * 
	 * This function is also used to calculate the Mi Matrix for the constant features.
	 */
	void calcMi(double[] mi, int overallPosition, int current, int stop, boolean doExp) {
		// Features will always be in the order of the potentials. We loop through the potentials, grabbing the features
		// for each.
		byte cachedPotential = -1;
		double cachedVal = Double.NaN;
		if (current < stop) {
			cachedPotential = potentialIx[current];
			cachedVal = val[current] * lambda[id[current]];
			//log.info(String.format("Pos: %d - pot %d, current: %d stop: %d id: %d", overallPosition, cachedPotential, current, stop, id[current]));
			++current;
		}
		double nodeVal = Double.NaN;
		int invalidIndex = overallPosition * nPotentials;
		for (short potential : orderedPotentials) {
			boolean invalid = overallPosition != -1 && invalidTransitions[invalidIndex + potential];
			double features = invalid ? Double.NEGATIVE_INFINITY : 0.0;
			// Add up all features for this potential.
			while (cachedPotential == potential) {
				features += cachedVal;
				// Assert.a(!Double.isNaN(features));
				if (current < stop) {
					cachedVal = val[current] * lambda[id[current]];
					cachedPotential = potentialIx[current];
					//log.info(String.format("Pos: %d - pot %d, current: %d stop: %d id: %d", overallPosition, cachedPotential, current, stop, id[current]));
					++current;
				} else {
					break;
				}
			}
			if (potential < nStates) {
				nodeVal = features;
			} else {
				int transition = potential - nStates;
				double val1 = features + nodeVal + constMi[transition];
				if (doExp) {
					val1 = exp(val1);
				}
				mi[transition] = val1;
			}
		}
		// verify that we actually checked all the cache entries.
		Assert.a(current == stop, "Pos: ", overallPosition, " Expected ", stop, " features only found ", current, " Pot: ", cachedPotential, " Val: ", cachedVal);
	}

	abstract class LengthFeatureProcessor {
		abstract void potential(short lookback, byte toNode, int transition, float value, int nodeStart, int nodeStop, int edgeStart, int edgeStop);

		/**
		 * Used to iterate through the cache of length based features and calculate weighed potential values for
		 * different lookback sizes and potentials. Used for updates in both the alpha and beta vectors.
		 * 
		 * Because the cache is so condensed, this function is very tricky.  The whole idea is that we rely on the cache being very well ordered
		 */
		void lengthCache(int overallPosition) {
			// When going through the length cache, you have to monitor for transitions between different lookbacks and
			// potentials
			short lastLookback = -1;
			byte lastPotential = -1;
			byte toNode = -1;
			float nodeValue = 0.0f;
			float potentialValue = 0.0f;
			int lengthCacheStart = lengthStarts[overallPosition];
			int lengthCacheStop = lengthStarts[overallPosition + 1];
			int nodeStart = -1;
			int nodeStop = -1;
			int edgeStart = -1;
			while (lengthCacheStart < lengthCacheStop) {
				// Get the next cache entry
				short lookback = lookbacks[lengthCacheStart];
				short index = lengthIndexes[lengthCacheStart];
				byte lengthPotential = lengthPotentials[lengthCacheStart];
				// See if we changed lookback or potential
				if (lookback != lastLookback || lengthPotential != lastPotential) {
					if(lengthPotential < nStates) {
						// Our current potential is a node and it differs from our previous potential.
						if(lastPotential != -1) {
							if(lastPotential < nStates) {
								// This is a node followed by a node, the first node must have been a valid starting node.
								//log.debug("Node followed by node");
								potential(lastLookback, lastPotential, -1, potentialValue, nodeStart, lengthCacheStart, -1, -1);
							}
							else {
								// This is an edge followed by a node.  Send out the previous edge.
								// It was an edge. update the appropriate previous matrices
								//log.debug("Edge followed by node");
								potential(lastLookback, toNode, lastPotential - nStates, potentialValue, nodeStart, nodeStop, edgeStart, lengthCacheStart);
							}
						}
						
						// Store the relevant stuff for this node.
						toNode = lengthPotential;
						nodeValue = 0.0f;
						potentialValue = 0.0f;
						nodeStart = lengthCacheStart;
						nodeStop = -1;
					}
					else {
						// This is an edge.
						// Deal with the last potential
						byte node = (byte) transitionTo[lengthPotential - nStates];
						if(lastPotential == -1) {
							toNode = node;
						}
						else {
							if (lastPotential < nStates) {
								// Last completed potential was a node.  check if it is related to this edge or not.
								if (lookback == lastLookback && node == toNode) {
									//log.debug("Node followed by related edge");
									// Node features are retrieved from the cache and updated for each edge they are a part of.
									// This was the end of a node's features.  Mark it.
									nodeStop = lengthCacheStart;
									nodeValue = potentialValue;
								}
								else {
									// The last potential was a node that had no edges.  Must have been a starting node
									//log.debug("node followed by unrelated edge");
									potential(lastLookback, toNode, -1, potentialValue, nodeStart, lengthCacheStart, -1, -1);

									// Reset our node value.
									nodeValue = 0.0f;
									nodeStart = -1;
									nodeStop = -1;
									toNode = node;
								}
							}
							else {
								//log.debug("Edge/edge.");
								// Last completed node was an edge.
								// Update the appropriate previous matrices
								int lastTrans = lastPotential - nStates;
								byte lastNode = (byte) transitionTo[lastPotential - nStates];
								potential(lastLookback, lastNode, lastTrans, potentialValue, nodeStart, nodeStop, edgeStart, lengthCacheStart);
								
								// Check if it was an edge for the same node.
								// Reset our node value if we changed nodes.
								if(lastNode != toNode || lookback != lastLookback) {
									nodeValue = 0.0f;
									nodeStart = -1;
									nodeStop = -1;
									toNode = node;
								}
							}
						}
						
						// Now set up the current edge
						edgeStart = lengthCacheStart;
						potentialValue = nodeValue;
					}

					// Update our last pointers to point to the current entry
					lastLookback = lookback;
					lastPotential = lengthPotential;
				}
				// Update the current potential total with a weighted feature sum
				if (index != -1) {
					potentialValue += lengthVals[lengthCacheStart] * lambda[index];
				}
				++lengthCacheStart;
			}
			if (lastPotential >= nStates) {
				//log.debug("Ending with edge");
				potential(lastLookback, toNode, lastPotential - nStates, potentialValue, nodeStart, nodeStop, edgeStart, lengthCacheStart);
			}
			else if(lastPotential != -1) {
				//log.debug("Ending with node");
				// If the last potential was a node, then this is a starting state.
				potential(lastLookback, lastPotential, -1, potentialValue, nodeStart, lengthCacheStart, -1, -1);
			}
		}
	}

	/**
	 * Given a vector with an existing normalization factor, convert it to a new normalization factor by scaling the
	 * entries.
	 */
	void renormalize(double[] vec, int currentNorm, int newNorm) {
		// Instead of dividing by the different (new-current), we reverse the subtraction to negate the exponent and
		// then multiply.
		double factor = exp(NORM_FACTOR * (currentNorm - newNorm));
		int len = vec.length;
		for (int i = 0; i < len; ++i) {
			vec[i] *= factor;
		}
	}

	/** Given a vector, computes a normalization factor for the entries and scales them according to that factor. */
	int normalize(double[] vec) {
		double sum = 0.0;
		int len = vec.length;
		for (int i = 0; i < len; ++i) {
			sum += vec[i];
		}
		if(sum == 0.0) {
			return 0;
		}
		Assert.a(!Double.isNaN(sum));
		if (sum > NORM_MIN && sum < NORM_MAX) {
			// No normalization required, our vector is in range.
			return 0;
		}
		//log.info("performing normalization");
		double val1 = log(sum);
		int norm = (int) val1 / NORM_FACTOR;
		val1 = exp(NORM_FACTOR * norm);
		for (int i = 0; i < len; ++i) {
			vec[i] /= val1;
		}
		return norm;
	}
}
