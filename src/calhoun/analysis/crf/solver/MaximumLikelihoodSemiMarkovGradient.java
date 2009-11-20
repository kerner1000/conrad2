package calhoun.analysis.crf.solver;

import java.io.BufferedWriter;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.CRFObjectiveFunctionGradient;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.analysis.crf.solver.CacheProcessor.FeatureEvaluation;
import calhoun.analysis.crf.solver.CacheProcessor.LengthFeatureEvaluation;
import calhoun.analysis.crf.solver.CacheProcessor.SolverSetup;
import calhoun.util.Assert;
import calhoun.util.ColtUtil;
import calhoun.util.FileUtil;


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
//public class MaximumLikelihoodSemiMarkovGradient extends CleanMaximumLikelihoodSemiMarkovGradient {
//}
public class MaximumLikelihoodSemiMarkovGradient implements CRFObjectiveFunctionGradient {
	private static final Log log = LogFactory.getLog(MaximumLikelihoodSemiMarkovGradient.class);
	private static final boolean debug = log.isDebugEnabled();
	private static final double ASSERTION_TOLERANCE = 0.0001;
	
	private static final int NORM_FACTOR = 50;
	private static final double NORM_MIN = Math.exp(-NORM_FACTOR);
	private static final double NORM_MAX = Math.exp(NORM_FACTOR);

	String alphaFile = null;
	String alphaLengthFile = null;
	String betaLengthFile = null;
	String expectFile = null;
	String expectLengthFile = null;
	String nodeMarginalFile = null;
	BufferedWriter alphaWriter = null;
	BufferedWriter alphaLengthWriter = null;
	BufferedWriter betaLengthWriter = null;
	BufferedWriter expectWriter = null;
	BufferedWriter expectLengthWriter = null;
	BufferedWriter nodeMarginalWriter = null;
	
	SolverSetup modelInfo;
	CacheProcessor cacheProcessor;
	FeatureEvaluation[] evals;
	LengthFeatureEvaluation[][] lengthEvals;
	boolean[] invalidTransitions;
	
	// / Cache feature information
	// / Cached value of the Mi matrix for all of the features present at every position
	// / Mi is stored as a sparse matrix
	short maxLookback;
	CacheProcessor.StatePotentials[] statesWithLookback;
	CacheProcessor.StatePotentials[] statesWithoutLookback;
	List<TrainingSequence> data;
	int iter = 0;
	double[][] alphas;
	int[] alphaNorms;
	double[] starterAlpha;

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

		// Initialize betas (for use later, in the gradient computation)
		alphas = new double[modelInfo.longestSeq][modelInfo.nStates];
		alphaNorms = new int[modelInfo.longestSeq];
		expects = new double[modelInfo.nFeatures];

		LookbackBuffer[] bufferContents = new LookbackBuffer[maxLookback+3];
		for(int i = 0; i<maxLookback+3; ++i) {
			bufferContents[i] = new LookbackBuffer();
		}
		lookbackBuffer = new RecyclingBuffer<LookbackBuffer>(bufferContents);
		nextBuffer = new LookbackBuffer();
		
		alphaProcessor = new AlphaLengthFeatureProcessor(); 
		betaProcessor = new BetaLengthFeatureProcessor(); 

		starterAlpha = new double[modelInfo.nStates];
	}

	public double apply(double[] param, double[] grad) {
		log.debug(String.format("Beginning It: %d Weights: %s", iter, ColtUtil.format(param)));
		alphaWriter = FileUtil.safeOpen(alphaFile);
		alphaLengthWriter = FileUtil.safeOpen(alphaLengthFile);
		betaLengthWriter = FileUtil.safeOpen(betaLengthFile);
		expectWriter = FileUtil.safeOpen(expectFile);
		expectLengthWriter = FileUtil.safeOpen(expectLengthFile);
		nodeMarginalWriter = FileUtil.safeOpen(nodeMarginalFile);
		lambda = param;
		Arrays.fill(grad, 0);
		double totalZ = 0.0;
		double result = 0.0;

		try {
			// Iterate through sequences
			Arrays.fill(expects, 0);
			for (int i = 0; i < modelInfo.nSeqs; ++i) {
				int len = modelInfo.seqOffsets[i + 1] - modelInfo.seqOffsets[i];
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
				//log.debug("Seq: "+i+" Z: "+printNorm(1/zInv, zNorm));
	
				// Work backwards, computing betas and expectations.
				betaProcessor.computeBetasAndExpectations(i, len);
	
				if(log.isDebugEnabled()) {
					double[][] seqFeatureSums = cacheProcessor.getSequenceFeatureSums();
					if(seqFeatureSums != null) {
						double seqResult = 0.0;
						for (int j = 0; j < modelInfo.nFeatures; ++j) {
							seqResult += seqFeatureSums[i][j] * param[j];
						}
						log.debug(String.format("Seq: %d L: %g LL: %f Training path: %f Z: %f", i, exp(seqResult-logZ), seqResult-logZ, seqResult, logZ));
						Assert.a(exp(seqResult-logZ) < 1.0);
					}
				}
				
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
				if(log.isDebugEnabled()) {
					log.debug(String.format("It: %d L=%e, LL=%f, norm(grad): %f Sums: %s Expects: %s Weights: %s Grad: %s", iter, exp(result), result,
							ColtUtil.norm(grad), ColtUtil.format(featureSums), ColtUtil.format(expects), ColtUtil.format(param), ColtUtil.format(grad)));
				}
				else {
					log.info(String.format("It: %d LL=%f, norm(grad): %f", iter, exp(result), result, ColtUtil.norm(grad)));
				}
			}
			Assert.a(exp(result) <= 1.0, "Likelihood is greater than 1.");
			result = result/modelInfo.totalPositions; 
			for(int i=0; i<grad.length; ++i) {
				grad[i] = grad[i]/modelInfo.totalPositions;
			}
			iter += 1;
		}
		finally {
			FileUtil.safeClose(alphaWriter);
			FileUtil.safeClose(alphaLengthWriter);
			FileUtil.safeClose(betaLengthWriter);
			FileUtil.safeClose(expectWriter);
			FileUtil.safeClose(expectLengthWriter);
			FileUtil.safeClose(nodeMarginalWriter);
		}
			
		return result;
	}

	public void clean() {
		// Clean up as much as possible
		modelInfo = null;
		cacheProcessor = null;
		evals = null;
		lengthEvals = null;
		invalidTransitions = null;
		statesWithLookback = null;
		statesWithoutLookback = null;
		data = null;
		alphas = null;
		alphaNorms = null;
		starterAlpha = null;
		lookbackBuffer = null;
		nextBuffer = null;
		expects = null;
		alphaProcessor = null; 
		betaProcessor = null; 
	}

	private final class AlphaLengthFeatureProcessor {
		int seqOffset;
		int pos;
		double[] alpha;
		int alphaNorm;
		double[] stableState;
		
		/**
		 * In the forward pass we compute alpha values and expections. This is simpler than the backwards pass because
		 * the cache is set up for us so that we can always look at one position at a time. We have to cache previous
		 * values but we never have to look ahead.
		 */
		final void computeAlpha(final int seqNum, final int len) {
			// Result the alpha norms
			Arrays.fill(alphaNorms, Integer.MIN_VALUE);

			Arrays.fill(starterAlpha, 0.0);
			double[] prevAlpha = null;
			
			seqOffset = modelInfo.seqOffsets[seqNum];

			for(pos = 0; pos < len; ++pos) {
				prevAlpha = alpha;
				alpha = alphas[pos];
				Arrays.fill(alpha, 0.0);
				if (pos == 0) {
					alphaNorm = 0;
					calcStartAlpha(alpha, seqNum);

					// Put an empty entry in the lookback so the first base has 0's initialized.
					Arrays.fill(nextBuffer.stableState, 0.0);
				} else {
					cacheMi(seqNum, nextBuffer.mi, stableState, nextBuffer.stableState, pos);
					regularAlphaUpdate(pos, nextBuffer.mi, prevAlpha, alpha);
				}
				
				// Add the lookback into the array
				stableState = nextBuffer.stableState;
				nextBuffer = lookbackBuffer.addFirst(nextBuffer);
					
				// Now we need to loop through for the length dependent cache
				lengthAlpha(seqNum, pos);

				int norm = normalize(alpha);
				/*if(norm != 0) {
					log.info("Pos: "+pos+" Renormalized alpha by "+norm+" to "+(alphaNorm+norm)+" : "+ColtUtil.format(alpha));
				}*/
				alphaNorm += norm; 
				alphaNorms[pos] = alphaNorm;
//				if(alphaNorm > 0 && pos > 0 && alphaNorms[pos-1]<=0)
//					log.info("Norm prob at pos: "+pos);
			}
		}

		/**
		 * Updates the alpha vector for non-length dependent states. We don't have to worry about normalization here
		 * because the regular alpha update is done before the length dependent, so these are the first values that will
		 * be set.
		 */
		private final void regularAlphaUpdate(final int pos, final double[] mi, final double[] lastAlpha, final double[] newAlpha) {
			double nodeVal = 0.0;
			int lastState = -1;
			boolean lengthNode = false;
			for (short potential : modelInfo.orderedPotentials) {
				if (potential < modelInfo.nStates) {
					if (lastState != -1) {
						newAlpha[lastState] = nodeVal;
					}
					lastState = potential;
					nodeVal = 0.0;
					lengthNode = modelInfo.maxStateLengths[potential] > 1;
				} else {
					if (!lengthNode) {
						int trans = potential - modelInfo.nStates;
						double transVal = mi[trans];
						if(!Double.isInfinite(transVal)) {
							int from = modelInfo.transitionFrom[trans];
							if(alphaWriter != null)
								FileUtil.safeWrite(alphaWriter, String.format("alpha[%d][%d] = %s = %s + alpha[%d][%d] %s * %s exp(%f)\n", pos, lastState, printNorm(nodeVal + lastAlpha[from] * exp(mi[trans]), alphaNorm), printNorm(nodeVal, alphaNorm), pos-1, from, printNorm(lastAlpha[from], alphaNorm), printNorm(exp(mi[trans]), 0), mi[trans]));
							nodeVal += lastAlpha[from] * exp(transVal);
						}
					}
				}
			}
			newAlpha[lastState] = nodeVal;
		}

		/** Updates an alpha entry with a weighted sum of features values for a given potential */
		private final void lengthAlpha(final int seqNum, final int pos) {
			cacheProcessor.evaluateSegmentsEndingAt(seqNum, pos);
			/*
			 * Updates an existing alpha by adding in: potentialValue - The value of any length-dependent features for
			 * this node f(y, i, d) and edge f(y', y, i, d) stableValue - The value of the non-length dependent node
			 * features summed across the length of this segment mis - The value of the non-length dependent transition
			 * from the previous node to this one f(y', y, i-d)
			 */
			int nSemiMarkovStates = modelInfo.statesWithLookback.length;
			for(int i=0; i<nSemiMarkovStates; ++i) {
				LengthFeatureEvaluation[] lookbacksForState = lengthEvals[i];
				CacheProcessor.StatePotentials statePotentials = modelInfo.statesWithLookback[i];
				byte toNode = statePotentials.state;
				
				int lbIndex=0;
				LengthFeatureEvaluation lengthEval = lookbacksForState[lbIndex];
				int lookback = lengthEval.lookback;
				while(lookback != -1) {
					//log.info("Pos: "+pos+"\t State: "+modelInfo.statesWithLookback[i].state+"\t Lookback: "+lookback);
					int prevPos = pos - lookback - 1;
					// For speed I hand inline RecyclingBuffer.get
					LookbackBuffer buffer = lookbackBuffer.array[(lookbackBuffer.currentStart+lookback)%lookbackBuffer.length];

					// Handle evaluation of the node potentials
					FeatureEvaluation nodeEvals = lengthEval.nodeEval;
					short[] indices = nodeEvals.index;
					float[] vals = nodeEvals.value;
					int ix = 0;
					short index = indices[ix];
					double stableValue = stableState[toNode] - buffer.stableState[toNode];
					double nodePotential = stableValue;
					while(index >= 0) {
						nodePotential += vals[ix] * lambda[index];
						index = indices[++ix];
					}
					if(debug)
						Assert.a(index != Short.MIN_VALUE, "Node lengths should only be returned in the cache if they are valid");

					if(prevPos < 0) {
						double nodeVal = nodePotential + starterAlpha[toNode];
						// If this is the first segment, then we don't worry about edges and handle the node directly.
						int norm = ((int) nodeVal) / NORM_FACTOR;
						nodeVal -= norm * NORM_FACTOR;
			
						if (norm > alphaNorm) {
							renormalize(alpha, alphaNorm, norm);
							//log.info("Renormalized alpha: "+ColtUtil.format(alpha));
							alphaNorm = norm;
						} else if (norm < alphaNorm) {
							nodeVal += NORM_FACTOR * (norm - alphaNorm);
						}
						if(alphaLengthWriter != null) {
							FileUtil.safeWrite(alphaLengthWriter, String.format("seq: %d alpha[%d][%d] = %s = %s + %s (Pot: %f Starter: %f)\n", seqNum, pos, toNode, printNorm(alpha[toNode] + exp(nodeVal), alphaNorm), printNorm(alpha[toNode], alphaNorm), printNorm(exp(nodeVal), alphaNorm), nodePotential, starterAlpha[toNode])); 
						}
							/*if((pos == 499 && toNode == 0))
							log.info(String.format("alpha[%d][%d] = %s = %s + %s (Pot: %f Starter: %f)", pos, toNode, printNorm(alpha[toNode] + exp(nodeVal), alphaNorm), printNorm(alpha[toNode], alphaNorm), printNorm(exp(nodeVal), alphaNorm), nodePotential, starterAlpha[toNode]));
							*/
						alpha[toNode] += exp(nodeVal);
					}
					else {
						// If this is not the first segment, we need to deal with edges coming into this segment
						FeatureEvaluation[] edgeEvals = lengthEval.edgeEvals;
						int nEdges = statePotentials.potentials.length;
						for(int edgeIx=0; edgeIx < nEdges; ++edgeIx) {
							int potential = statePotentials.potentials[edgeIx];
							int trans = potential - modelInfo.nStates;
							int fromNode = modelInfo.transitionFrom[trans];
							// Skip semi-Markov self transitions
							if(fromNode == toNode)
								continue;

							double edgeVal = 0.0;

							if(edgeEvals == null) {
								// If the cache processor does not have edge evaluations
								// Just check if this transition is legal based on the invalid transitions matrix
								int invalidIndex = (seqOffset+prevPos+1)*modelInfo.nPotentials;
								if(invalidTransitions[invalidIndex + potential]) {
									//log.info("Illegal transition: "+fromNode+"-"+toNode+" at pos: "+prevPos);
									continue;
								}
							}
							else {
								// If the cache processor does have edge evaluations, then ignore the illegal transitions matrix
								// and update the expval using the edge evaluations
								FeatureEvaluation potEvals = edgeEvals[edgeIx];
								indices = potEvals.index;
								vals = potEvals.value;
								ix = 0;
								index = indices[i];
								while(index >= 0) {
									edgeVal += vals[ix] * lambda[index];
									index = indices[++ix];
								}
								if(index == Short.MIN_VALUE) {
									continue;
								}
							}
							
							// Renormalize and update the exp value.
							double expVal = edgeVal + buffer.mi[trans] + nodePotential;
							int expNorm = ((int) expVal)/NORM_FACTOR;
							expVal -= expNorm*NORM_FACTOR;
							
							int prevNorm = alphaNorms[prevPos];
							int updateNorm = expNorm + prevNorm;
							if(updateNorm > alphaNorm) {
								// Our updated value is larger than the existing alpha value, renormalize that alpha vector.
								renormalize(alpha, alphaNorm, updateNorm);
								//log.info("Renormalized alpha from "+alphaNorm+" to "+updateNorm +" : "+ColtUtil.format(alpha));
								alphaNorm = updateNorm;
							}
							else if(alphaNorm > updateNorm) {
								// Renormalize the expVal 
								int expShift = alphaNorm - updateNorm;
								//log.info(String.format("Renormalize feature by %d from %d to %d",expShift, expNorm, expNorm+expShift));
								expNorm += expShift;
								expVal -= expShift*NORM_FACTOR;
							}
							
							double prevAlpha = alphas[prevPos][fromNode];
							double update = exp(expVal) * prevAlpha;
							if(alphaLengthWriter != null) {
								FileUtil.safeWrite(alphaLengthWriter, String.format("seq: %d alpha[%d][%d] = %s = %s + %s (alpha[%d][%d]) * %s exp(EdgeLength: %f NodeLength: %f Edge: %f Node: %f )\n", 
										seqNum, pos, toNode, printNorm(alpha[toNode] + update, alphaNorm), printNorm(alpha[toNode], alphaNorm), printNorm(prevAlpha, alphaNorms[prevPos]), 
										prevPos, modelInfo.transitionFrom[trans], printNorm(exp(expVal), expNorm), edgeVal, nodePotential - stableValue, buffer.mi[trans], stableValue));
							}
										
							alpha[toNode] += update;
							// Expensive assertion that can catch some normalization problems
							if(debug)
								Assert.a(expNorm + prevNorm == alphaNorm, "Norm problem.  Exp: ", expNorm, " Prev alpha: ", prevNorm, " Alpha: ", alphaNorm);
						}
					}
					
					++lbIndex;
					lengthEval = lookbacksForState[lbIndex];
					lookback = lengthEval.lookback;
				}
			}
		}

		/**
		 * A specialized version of calcMi for the first position in a sequence. Has the special property that constant
		 * edge features are not included.
		 */
		void calcStartAlpha(double[] currentAlpha, int seq) {
			cacheProcessor.evaluatePosition(seq, 0);
			int invalidIndex = seqOffset*modelInfo.nPotentials;
			for(short potential : modelInfo.orderedPotentials) {
				if(potential < modelInfo.nStates) {
					boolean invalid = invalidTransitions[invalidIndex + potential];
					double features = invalid ? Double.NEGATIVE_INFINITY : 0.0;

					// Add up all features for this potential.
					FeatureEvaluation potEvals = evals[potential];
					short[] indices = potEvals.index;
					float[] vals = potEvals.value;
					int i = 0;
					short index = indices[i];
					while(index != -1) {
						if(index == Short.MIN_VALUE) {
							features = Double.NEGATIVE_INFINITY;
							break;
						}
						features += vals[i]*lambda[index];
						index = indices[++i];
					}
					if(modelInfo.maxStateLengths[potential]> 1) {
						starterAlpha[potential] = features;
					}
					else {
						currentAlpha[potential] = exp(features);
					}
				}
			}		
		}
	}

	/** Updates a beta entry with a potential from a length dependent feature */
	class BetaLengthFeatureProcessor {
		int seqOffset;
		int lengthPos;
		double[] lengthStable;
		int miPos;
		double[] stableState;

		double[] beta;
		int betaNorm;

		double prob;
		
		LookbackBuffer posLookback; 
		LookbackBuffer prevLookback; 
		
		// This vector tracks the marginal probability of each node, accounting for length dependence
		double[] nodeProb = new double[modelInfo.nStates];
		double[] newNodeProb = new double[modelInfo.nStates];

		// This vector tracks the marginal probability of each edge, accounting for length dependence
		double[] edgeProb = new double[modelInfo.nTransitions];
		
		/**
		 * Computing the betas is a bit tricky. The problem is that our cache associates each explicit
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
		void computeBetasAndExpectations(int seqNum, int len) {
			// Start by initializing an array of unexponentiated mi and stable state matrices for all of the possible
			// lookback positions from the last base
			// Also initialize the array of previous beta vectors
			seqOffset = modelInfo.seqOffsets[seqNum];

			// lastInitPos holds the position of the leftmost position which gets evaluated at the start of the run.
			int lastInitPos = len - 2 - maxLookback;

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
						nextBuffer.clear();
						cacheMi(seqNum, nextBuffer.mi, stableState, nextBuffer.stableState, miPos+1);
					}
					nextBuffer.pos = miPos;
					stableState = nextBuffer.stableState;
					nextBuffer = lookbackBuffer.addFirst(nextBuffer);
					nextBuffer.clear();
					--miPos;
				}

				// At this point, miPos contains the leftmost position for which markov features have been evaluated.
				// The lookback buffer at this position has a stableState vector for each length dependent state and
				// an mi matrix computed for that state.
				
				// Retrieve the lookback information for the current position.  The lookback buffer runs left to right, so the 
				// current position is not at the beginning.
				posLookback = lookbackBuffer.get(pos - miPos -1);

				// Check that posLookback is correct
				if(debug)
					Assert.a(posLookback.pos == pos, "Wrong lookback buffer: was ", posLookback.pos, " should be ", pos);

				// Update betas (if necessary)
				if(prevLookback == null) {
					// For the last position in the sequence (first iteration in the for loop)
					if(debug)
						Assert.a(pos == len-1);

					// Betas for the last position are all 1.  exp(0)
					Arrays.fill(posLookback.beta, 1);
					posLookback.betaNorm = 0;

					/* Initialize node marginals.  Since all segments end at the last position, we can compute
					this by multiplying the alpha and beta vectors and dividing by Z.  The beta vector is all 1's though
					so we just divide alpha by z */ 
					double nodeNorm = exp((alphaNorms[pos] - zNorm) * NORM_FACTOR) * zInv;
					for(int i=0; i<modelInfo.nStates; ++i) {
						nodeProb[i] = nodeNorm * alphas[pos][i];
					}
					log.debug("Node marginals at seq "+seqNum+" last position ("+pos+"): "+ColtUtil.format(nodeProb));
				}
				else {
					// For all but the last position, we update beta first for the markov states and then the semi-markov states.
					if(debug)
						Assert.a(prevLookback.pos == pos+1);
					posLookback.betaNorm = regularBetaUpdate(pos+1, posLookback.beta, posLookback.betaNorm, prevLookback.beta, prevLookback.betaNorm, prevLookback.transitionProb, posLookback.mi);
					//log.info(String.format("Pos: %d, Mipos: %d LbPos: %d", pos, miPos, posLookback.pos));

					// At this point, the nodeProb and edgeProb vectors contain the marginal probabilities of markov nodes and edges into them.
					// The beta is fully updated for the markov states.
					// Semi markov has not been done.
					
					// Now we need to loop through for the length dependent cache
					beta = prevLookback.beta;
					/*for(int ix=0; ix<beta.length; ++ix) {
						Assert.a(!Double.isNaN(beta[ix]), prevLookback.pos+ " "+ ix);
					}*/
					betaNorm = prevLookback.betaNorm;
					lengthStable = prevLookback.stableState;
					lengthPos = pos + 1;
					lengthBeta(seqNum, lengthPos);

					System.arraycopy(nodeProb, 0, newNodeProb, 0, modelInfo.nStates);
					
					// Now calculate edge marginals for staying in the same state
					// Take the node probabilities, and subtract off each transition
					for(CacheProcessor.StatePotentials lb: statesWithLookback) {
						int state = lb.state;
						int index = modelInfo.selfTransitions[state];
						double transProb = nodeProb[state];
						for(byte pot : lb.potentials) {
							double lbTrans = posLookback.transitionProb[pot - modelInfo.nStates];
							transProb -= lbTrans;
							newNodeProb[state] -= lbTrans;
						}
						Assert.a(posLookback.transitionProb[index] == 0.0);
						//log.debug("Self-trans marginal for ("+state+"): "+transProb );
						posLookback.transitionProb[index] = transProb;
					}
				}
				
				// As a check, we verify that the node marginals sum to one for each position.
				double sum =0.0;
				for(double x : nodeProb) {
					if(x > 1.0+ASSERTION_TOLERANCE || x < -ASSERTION_TOLERANCE)
						Assert.a(false, "Iter ",iter," Pos: "+pos+" Node marginals not valid "+x);
					sum += x;
				}
				if(Math.abs(1.0-sum) > ASSERTION_TOLERANCE) {
					Assert.a(false, "Iter ",iter," Pos: "+pos+" Node marginals sum to "+sum+" not 1: ", ColtUtil.format(nodeProb), " at ", seqNum," ",pos);
				}

				// Verify that edge marginals sum to 1.
				if(prevLookback != null) {
					if(debug) {
						sum =0.0;
						for(double x : edgeProb) {
							if(x > 1+ASSERTION_TOLERANCE || x < -ASSERTION_TOLERANCE)
								Assert.a(false, "Iter ",iter," Pos: "+pos+" Edge marginal not valid "+x);
							sum += x;
						}
						for(double x : posLookback.transitionProb) {
							if(x > 1+ASSERTION_TOLERANCE || x < -ASSERTION_TOLERANCE)
								Assert.a(false, "Iter ",iter," Pos: "+pos+" Self-trans marginal not valid "+x);
							sum += x;
						}
						/*
						for(int i=0; i<modelInfo.nStates; ++i) {
							log.debug("Nod "+i+" "+nodeProb[i]);
						}*/
						if(Math.abs(1.0-sum) > 0.001) {
							/*for(int i = 0; i < posLookback.transitionProb.length; ++i) {
								if(posLookback.transitionProb[i] != 0)
									log.debug("Seg "+modelInfo.transitionFrom[i]+"-"+modelInfo.transitionTo[i]+" "+posLookback.transitionProb[i]);
							}
							for(int i = 0; i < edgeProb.length; ++i) {
								if(edgeProb[i] != 0)
									log.debug("Edg "+modelInfo.transitionFrom[i]+"-"+modelInfo.transitionTo[i]+" "+edgeProb[i]);
							}*/
							Assert.a(false, "Edge marginals don't sum to 1.  Sum to: ", sum, " - ", ColtUtil.format(edgeProb), ColtUtil.format(posLookback.transitionProb));
						}
					}
					
					// At this point the previous beta values are all updated.
					// Update expectations for the transitions we just calculated.
					updateExpectations(seqNum, pos+1, posLookback.transitionProb);
				}

				// Update the node probabilities to remove the length transitions
				double[] temp = nodeProb;
				nodeProb = newNodeProb;
				newNodeProb = temp;
				
				/*if (debug) {
					if ((seqOffset == 0) && (pos < 2 || pos >= len - 2)) {
						log.debug(String.format("Pos: %d expects: %s alphas: %s (norm %d) betas: %s (norm %d) MiPos: %d", pos, ColtUtil.format(expects), ColtUtil
								.format(alphas[pos]), alphaNorms[pos], ColtUtil.format(posLookback.beta), posLookback.betaNorm, miPos+1));
					}
				}*/

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
			lengthBeta(seqNum, lengthPos);
			//betaProcessor.lengthCache(seqStart);
			updateExpectations(seqNum, 0, posLookback.transitionProb);
		}

		/**
		 * Does the update of the beta values for all of the regular, non-length dependent states. Also calculates all
		 * of the node and edge probabilities for the non-length nodes and the transitions into them.
		 */
		private int regularBetaUpdate(int pos, double[] newBeta, int newNorm, double[] oldBeta, int oldNorm, double[] transitionProb, double[] mi) {
			// Need to deal with possibly different normalization constants.
			int norm = newNorm;
			double normAdjust = 0.0;
			if (oldNorm > newNorm) {
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
				if(nodeMarginalWriter != null) {
					FileUtil.safeWrite(nodeMarginalWriter, String.format("NodeMarg[%d][%d] = %f = %f * %f * %f (aN: %d bN: %d zN: %d 1/z: %f)\n", pos, node, nodeAlpha[node] * betaVal * nodeNorm, nodeAlpha[node], betaVal, nodeNorm, alphaNorms[pos], oldNorm, zNorm, zInv));
				}
				nodeProb[node] = nodeAlpha[node] * betaVal * nodeNorm;
				
				// For regular states, we sum edge probabilities to get node probabilities.
				if(pos > 0) {
					for (short potential : potentials.potentials) {
						int trans = potential - modelInfo.nStates;
						double transVal = mi[trans];
						if(!Double.isInfinite(transVal)) {
							// Mi for beta is not exponentiated, so we do it here.
							double potentialValue = exp(mi[trans] + normAdjust);
							nodePotential += potentialValue;

							int from = modelInfo.transitionFrom[trans];
							/*if(pos > 2400) {
								log.info(String.format("Beta[%d][%d] = %f = %f + Beta[%d][%d] %f * %f (Mi: %f)", pos, from, newBeta[from] + potentialValue * betaVal, newBeta[from], pos+1, node, betaVal, potentialValue, mi[trans]));
							}*/
							newBeta[from] += potentialValue * betaVal;
							edgeProb[trans] = edgeAlpha[from] * potentialValue * betaVal * edgeNorm;
						}
						else {
							edgeProb[trans] = 0.0;
						}
					}
				}
			}

			// Now check to see if this new vector needs normalization again.
			//log.debug("pos "+pos);
			int ret = norm;
			if(newBeta != null) {
				try {
					ret += normalize(newBeta);
				}
				catch(RuntimeException ex) {
					log.warn("Normalization problem at " + pos + " "+ColtUtil.format(newBeta));
					throw ex;
				}
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
		private void lengthBeta(int seqNum, int pos) {
			cacheProcessor.evaluateSegmentsEndingAt(seqNum, pos);

			int nSemiMarkovStates = modelInfo.statesWithLookback.length;
			for(int i=0; i<nSemiMarkovStates; ++i) {
				LengthFeatureEvaluation[] lookbacksForState = lengthEvals[i];
				CacheProcessor.StatePotentials statePotentials = modelInfo.statesWithLookback[i];
				byte toNode = statePotentials.state;
				
				int lbArrayIndex=0;
				LengthFeatureEvaluation lengthEval = lookbacksForState[lbArrayIndex];
				int lookback = lengthEval.lookback;
				while(lookback != -1) {
					int prevPos = lengthPos - lookback - 1;
					int lbIndex = prevPos - miPos - 1;

					//log.info("Pos: "+pos+"\t State: "+modelInfo.statesWithLookback[i].state+"\t Lookback: "+lookback);
					LookbackBuffer segBegin = null;
					if(prevPos >= 0) {
						// For speed I hand inline RecyclingBuffer.get
						segBegin = lookbackBuffer.array[(lookbackBuffer.currentStart+lbIndex)%lookbackBuffer.length];
					}
					LookbackBuffer stableBuffer = lookbackBuffer.array[(lookbackBuffer.currentStart+lbIndex+1)%lookbackBuffer.length];

					// Handle evaluation of the node potentials
					double stableValue = stableBuffer.stableState[toNode] - lengthStable[toNode];
					double nodePotential = stableValue;
					
					FeatureEvaluation nodeEvals = lengthEval.nodeEval;
					short[] indices = nodeEvals.index;
					float[] vals = nodeEvals.value;
					int ix = 0;
					short index = indices[ix];
					while(index >= 0) {
						nodePotential += vals[ix] * lambda[index];
						index = indices[++ix];
					}
					if(debug)
						Assert.a(index != Short.MIN_VALUE, "Node lengths should only be returned in the cache if they are valid.  They can be invalid because a node is invalid or a self-transition edge is invalid.");

					if(prevPos < 0) {
						// If this is the first segment, then we don't worry about edges and handle the node directly.
						double expVal = nodePotential + starterAlpha[toNode];
						lengthBetaHandling(seqNum, prevPos, pos, expVal, -1, toNode, 1.0, 0, nodeEvals);
					}
					else {
						// If this is not the first segment, we need to deal with edges
						FeatureEvaluation[] edgeEvals = lengthEval.edgeEvals;
						int nEdges = statePotentials.potentials.length;
						for(int edgeIx=0; edgeIx < nEdges; ++edgeIx) {
							int potential = statePotentials.potentials[edgeIx];
							int trans = potential - modelInfo.nStates;
							int fromNode = modelInfo.transitionFrom[trans];
							// Skip semi-Markov self transitions
							if(fromNode == toNode)
								continue;

							double edgeVal = 0.0;

							if(edgeEvals == null) {
								// If the cache processor does not have edge evaluations
								// Just check if this transition is legal based on the invalid transitions matrix
								int invalidIndex = (seqOffset + prevPos+1)*modelInfo.nPotentials;
								if(invalidTransitions[invalidIndex + potential])
									continue;
							}
							else {
								// If the cache processor does have edge evaluations, then ignore the illegal transitions matrix
								// and update the expval using the edge evaluations
								FeatureEvaluation potEvals = edgeEvals[edgeIx];
								indices = potEvals.index;
								vals = potEvals.value;
								ix = 0;
								index = indices[i];
								while(index >= 0) {
									edgeVal += vals[ix] * lambda[index];
									index = indices[++ix];
								}
								if(index == Short.MIN_VALUE) {
									continue;
								}
							}

							//log.debug("Stab: "+ stableBuffer.stableState[toNode]+" - " + lengthStable[toNode]);
							if(debug) {
								Assert.a(prevLookback.pos == lengthPos, "Expected ",lengthPos, " was ",prevLookback.pos);
								Assert.a(segBegin.pos == (lengthPos-lookback-1), "Expected ",(lengthPos-lookback-1), " was ",segBegin.pos);
							}
							double expVal = edgeVal + segBegin.mi[trans] + nodePotential;
							double prevAlpha = alphas[prevPos][fromNode];
							int prevAlphaNorm = alphaNorms[prevPos];
							//log.debug("mi: "+buffer.mi[trans]+" s: "+nodePotential+" Base: "+(expVal - buffer.mi[trans] - nodePotential));
							int expNorm = lengthBetaHandling(seqNum, prevPos, pos, expVal, fromNode, toNode, prevAlpha, prevAlphaNorm, nodeEvals);
							expVal -= expNorm * NORM_FACTOR;

							// To determine the normalization value, we compare the existing beta value to the value we are about to add
							// Whichever is a larger number will dominate, and so we use that normalization value as the new value and ignore the
							// previous one.
							
							// Update the beta values
							int updateNorm = expNorm + betaNorm;
							if(updateNorm > segBegin.betaNorm) {
								//log.info(String.format("Renormalize beta from %d to %d", segBegin.betaNorm, updateNorm));
								renormalize(segBegin.beta, segBegin.betaNorm, updateNorm);
								segBegin.betaNorm = updateNorm;
							}
							else if(segBegin.betaNorm > updateNorm) {
								int expShift = updateNorm - segBegin.betaNorm;
								//log.info(String.format("Renormalize feature from %d to %d",expNorm, expNorm+expShift));
								expNorm += expShift;
								expVal += expShift*NORM_FACTOR;
							}
							double transPotential = exp(expVal);
							double update = transPotential * beta[toNode];
							segBegin.beta[fromNode] += update;

							// Update expectations
							if(edgeEvals != null) {
								FeatureEvaluation potEvals = edgeEvals[edgeIx];
								indices = potEvals.index;
								vals = potEvals.value;
								ix = 0;
								index = indices[i];
								while(index != -1) {
									if(expectLengthWriter != null)
										FileUtil.safeWrite(expectLengthWriter, String.format("Seq %d Pos %d-%d Expect #%d: %e = %e + Prob: %e * EdgeVal: %e\n", seqNum, prevPos, pos, index, expects[index]+prob*vals[i], expects[index], prob, vals[i]));
									expects[index] += prob * vals[ix];
									index = indices[++ix];
								}
							}
							
							// Updates the transition probabilities
							//log.info(String.format("EdgeMarg[%d][%d] = %f = %f + %f", prevPos, toNode, buffer.transitionProb[trans] + prob, buffer.transitionProb[trans], prob));
							segBegin.transitionProb[trans] += prob;

							if(betaLengthWriter != null) {
								FileUtil.safeWrite(betaLengthWriter, String.format(String.format("Beta[%d][%d] = %s = %s + %s beta[%d][%d] * %s exp(Edge: %f Node: %f Stable: %f Trans: %f)\n", 
										prevPos, fromNode, printNorm(segBegin.beta[fromNode], segBegin.betaNorm), printNorm(segBegin.beta[fromNode]-update, segBegin.betaNorm), printNorm(beta[toNode], betaNorm), lengthPos, toNode, printNorm(transPotential, expNorm), 
										edgeVal, nodePotential - stableValue, stableValue, segBegin.mi[trans]))); 
							}
						}
					}
					
					++lbArrayIndex;
					lengthEval = lookbacksForState[lbArrayIndex];
					lookback = lengthEval.lookback;
				}
			}
		}

		int lengthBetaHandling(int seqNum, int prevPos, int pos, double expVal, int fromNode, int toNode, double prevAlpha, int prevAlphaNorm, FeatureEvaluation nodeEvals) {
			int norm = ((int) expVal) / NORM_FACTOR;
			expVal -= norm * NORM_FACTOR;
			
			// In addition to updating the beta array, we need to calculate a probability for this segment so we can
			// correctly calculate feature expectations
			//log.info(String.format("PrevAlpha %f Beta: %f Exp: %f zInv: %f", prevAlpha, beta[toNode], expVal, zInv));
			double afterExp = exp(expVal + NORM_FACTOR * (prevAlphaNorm + norm + betaNorm - zNorm));
			// This if statement is here in case alpha or beta is 0, but the normalization is large.  This could cause the exp value to go to infinity and result in a NaN probability instead of 0.
			if(prevAlpha == 0 || beta[toNode] == 0)
				prob = 0.0;
			else
				prob = prevAlpha * beta[toNode] * afterExp * zInv;
			if(Double.isNaN(prob)) {
				log.info(String.format("NaN = Alpha: %s * Beta: %s * Seg: %s / Z: %s",printNorm(prevAlpha, prevAlphaNorm), printNorm(beta[toNode], betaNorm), printNorm(exp(expVal), norm), printNorm(1/zInv, zNorm)));
				Assert.a(false, String.format("Seq: %d Pos: %d-%d: Bad prob (NaN) = Alpha: %e * Beta[%d] %e * %e exp(%f Norm  a:%d n:%d b:%d z:%d) * %e",
						seqNum, prevPos, pos, prevAlpha, toNode, beta[toNode], afterExp, expVal, prevAlphaNorm, norm, betaNorm, zNorm, zInv));
			}
			
			// Now update expectations for all node features for this edge
			short[] indices = nodeEvals.index;
			float[] vals = nodeEvals.value;
			int i = 0;
			short index = indices[i];
			while(index != -1) {
				if(expectLengthWriter != null)
					FileUtil.safeWrite(expectLengthWriter, String.format("Seq %d Pos %d-%d Expect #%d: %e = %e + Prob: %e * NodeVal: %e\n", seqNum, prevPos, pos, index, expects[index]+prob*vals[i], expects[index], prob, vals[i]));
				if(prob != 0.0)
					expects[index] += prob * vals[i];
				index = indices[++i];
			}
			
			if(nodeMarginalWriter != null) {
				FileUtil.safeWrite(nodeMarginalWriter, String.format(
						"NodeMarg[%d][%d] = %f = %f + Alpha[%d][%d]: %s * Beta[%d][%d]: %s * seg: %s / Z: %s\n", 
						lengthPos, toNode, nodeProb[toNode]+prob, nodeProb[toNode], prevPos, fromNode, printNorm(prevAlpha, prevAlphaNorm), 
						pos, toNode, printNorm(beta[toNode], betaNorm), 
						printNorm(expVal, norm), printNorm(1/zInv, zNorm)));
			}
			nodeProb[toNode] += prob;
			
			return norm;
		}
		
		void updateExpectations(int seqNum, int pos, double[] transitionProb) {
			cacheProcessor.evaluatePosition(seqNum, pos);

			// First compute the expectations for the length-dependent states.
			int invalidIndex = (seqOffset+pos)* modelInfo.nPotentials;

			boolean lengthNode = false;
			for (short potential : modelInfo.orderedPotentials) {
				boolean invalid = invalidTransitions[invalidIndex + potential];
				
				double prob = Double.NaN;
				if(potential < modelInfo.nStates) {
					lengthNode = modelInfo.maxStateLengths[potential] > 1;
					prob = nodeProb[potential];
				}
				else {
					if(pos == 0)
						continue;
					prob = (lengthNode ? transitionProb : edgeProb)[potential-modelInfo.nStates]; 
					//prob = transitionProb[potential-nStates]; 
				}
				
				if (!invalid) {
					// Iterate through features for this potential.
					FeatureEvaluation potEvals = evals[potential];
					short[] indices = potEvals.index;
					float[] vals = potEvals.value;
					int i = 0;
					short index = indices[i];
					while(index != -1) {
						if(expectWriter != null)
							FileUtil.safeWrite(expectWriter, String.format("Seq %d Pos %d Expect #%d: %e = %e + Prob: %e * Val: %e\n", seqNum, pos, index, expects[index]+prob*vals[i], expects[index], prob, vals[i]));
						expects[index] += prob*vals[i];
						index = indices[++i];
					}
				}
			}
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
	private static final void renormalize(final double[] vec, final int currentNorm, final int newNorm) {
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
	private static final int normalize(final double[] vec) {
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

	private static final double exp(final double val) {
		return Math.exp(val);
	}

	private static final double log(final double val) {
		return Math.log(val);
	}

	/**
	 * This object holds information about previous positions during the computation of betas and expectations. This
	 * allows us to quickly access data about previous positions. These objects are kept in a recycling buffer that
	 * keeps one buffer for each possible lookback.
	 * 
	 * One tricky aspect of this is that the details change slightly between the forward and backwards pass.  On the forward
	 * pass, the lookback contains the information in the normal way.  In the backwards pass, stable states and transitions are 
	 * shifted back one base compared to the betas.
	 */
	final class LookbackBuffer {
		int pos;
		
		// The mi matrix for transitioning from pos-lookback-1 to pos-lookback
		double[] mi = new double[modelInfo.nPotentials];

		// The weighted sum of feature values for staying in this position from the beginning to pos-lookback
		double[] stableState = new double[modelInfo.nStates];

		// Initial values of the beta vector for somelength dependent states.
		double[] beta = new double[modelInfo.nStates];
		
		// Norm of prevBeta.
		int betaNorm;
		
		// Stores the probability of all segments begining at this position using this transition.
		double[] transitionProb = new double[modelInfo.nTransitions];

		/** mi and stableStates are cleared as new values are entered. This fixes the others */
		void clear()
		{
			pos = -1;
			Arrays.fill(beta, 0.0);
			betaNorm = Integer.MIN_VALUE;
			Arrays.fill(transitionProb, 0.0);
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
		return alphaLengthFile;
	}

	public void setAlphaLengthFile(String alphaLengthFile) {
		this.alphaLengthFile = alphaLengthFile;
	}

	public String getAlphaFile() {
		return alphaFile;
	}

	public void setAlphaFile(String alphaFile) {
		this.alphaFile = alphaFile;
	}

	public String getExpectFile() {
		return expectFile;
	}

	public void setExpectFile(String expectFile) {
		this.expectFile = expectFile;
	}

	public String getExpectLengthFile() {
		return expectLengthFile;
	}

	public void setExpectLengthFile(String expectLengthFile) {
		this.expectLengthFile = expectLengthFile;
	}

	public String getNodeMarginalFile() {
		return nodeMarginalFile;
	}

	public void setNodeMarginalFile(String nodeMarginalFile) {
		this.nodeMarginalFile = nodeMarginalFile;
	}

	public String getBetaLengthFile() {
		return betaLengthFile;
	}

	public void setBetaLengthFile(String betaLengthFile) {
		this.betaLengthFile = betaLengthFile;
	}

	public double[] getFeatureSums() {
		return this.featureSums.clone();
	}
}
