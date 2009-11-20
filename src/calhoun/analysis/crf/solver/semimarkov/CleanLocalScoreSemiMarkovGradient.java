package calhoun.analysis.crf.solver.semimarkov;

import java.util.Arrays;
import java.util.List;

import calhoun.analysis.crf.LocalPathSimilarityScore;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.analysis.crf.scoring.SimScoreMaxStateAgreement;
import calhoun.analysis.crf.solver.CacheProcessor;
import calhoun.analysis.crf.solver.LookbackBuffer;
import calhoun.analysis.crf.solver.CacheProcessor.FeatureEvaluation;
import calhoun.analysis.crf.solver.CacheProcessor.LengthFeatureEvaluation;
import calhoun.util.Assert;
import calhoun.util.ColtUtil;
import calhoun.util.FileUtil;

/** computes an objective function which is the expected value of a local path similarity score on a 
 * semi-Markov model.  Requires a {@link CacheProcessor} and a {@link LocalPathSimilarityScore} to be configured.<p>
 * <h2>Debugging output</h2>
 * To get a better understanding of what the objective function is doing, several different properties can be set that
 * cause the objective function to write out trace files showing its calculations during training.  Usually when turning
 * these options on, you should set <code>maxIters = 1</code> and <code>requireConvergence = false</code> in your optimizer
 * to do only a single training iteration, possibly setting the starts to some predetermined value.  Each of these
 * properties can be configured with a filename and each time {@link #apply} is called, the file will be overwritten with 
 * data from the current call.  The logging options are:
 *
 * <ul>
 * <li> <b><code>alphaFile</code></b> - computation of alpha values for Markov states, includes all nodes and edges.
 * <li> <b><code>alphaLengthFile</code></b> - computation of alpha values for semi-Markov states , includes all segments
 * <li> <b><code>betaLengthFile</code></b> - computation of beta values for semi-Markov states , includes all segments
 * <li> <b><code>expectFile</code></b> - computation of expected values for each Markov feature 
 * <li> <b><code>expectLengthFile</code></b> - computation of expected values for each semi-Markov feature  
 * <li> <b><code>nodeMarginalFile</code></b> - computation of marginal probability of each state at each position 
 * </ul>

 * <h4>Implementation Notes</h4>
 * The general normalization scheme works as follows. When updating alpha values in the forward pass we compute segments
 * of length 1 first and then work backwards.
 * <p>
 * Instead of always normalizing to 1 we discretize the normalization. We choose an arbitrary normalization factor w,
 * such as 50. The normalization factor at any position is then an integer v, and all entries at that position are
 * alpha[y]*e^(v*w).
 * <p>
 * The normalization can be computed at any position from 1) Elements of the alpha array are summed s 2) v = log(s)/w.
 * By integer division v will always be an appropriate normalizer. It may be positive or negative. 3) All elements of
 * the array are divided by e^(v*w)
 * 
 */
public class CleanLocalScoreSemiMarkovGradient extends CleanMaximumLikelihoodSemiMarkovGradient {
	LocalPathSimilarityScore score = new SimScoreMaxStateAgreement();

	// Score cache
	double[][] localScoreStableCache;
	double[][] localScoreTransitionCache;

	double[][] betas;
	int[] betaNorms;
	double[][] allEdgeProb;
	double[][] allNodeProb;

	double[][] scoreAlpha;
	double[][] scoreBeta;
	double[][] semiMarkovScoreAlpha;
	double[][] semiMarkovScoreBeta;

	boolean semiMarkov;
	
	@Override
	public void setTrainingData(ModelManager fm, List<? extends TrainingSequence<?>> data) {
		super.setTrainingData(fm, data);

		// Initialize arrays to hold betas and marginals
		betas = new double[modelInfo.nStates][modelInfo.longestSeq];
		betaNorms = new int[modelInfo.longestSeq];
		allEdgeProb = new double[modelInfo.nTransitions][modelInfo.longestSeq]; //[pos][yprev][y], not defined for pos=0
		allNodeProb = new double[modelInfo.nStates][modelInfo.longestSeq];
		scoreAlpha = new double[modelInfo.nStates][modelInfo.longestSeq];
		scoreBeta = new double[modelInfo.nStates][modelInfo.longestSeq];

		semiMarkov = modelInfo.maxLookback > 1;
		if(semiMarkov) {
			semiMarkovScoreAlpha = new double[modelInfo.nStates][modelInfo.longestSeq];
			semiMarkovScoreBeta = new double[modelInfo.nStates][modelInfo.longestSeq];
		}

		betaProcessor.setGlobalArrays(betas, betaNorms, allNodeProb, allEdgeProb);
		fillScoreCache(cacheProcessor.getData());
	}

	@Override
	public double apply(double[] param, double[] grad) {
		log.debug(String.format("Beginning It: %d Weights: %s", iter, ColtUtil.format(param)));
		lambda = param;
		double[] scoreFeatureProductExpectation = new double[grad.length];
		Arrays.fill(grad, 0);
		double result = 0.0;

		try {
			// Iterate through sequences
			logs.open();
			for (int i = 0; i < modelInfo.nSeqs; ++i) {
				Arrays.fill(expects, 0);

				int len = modelInfo.seqOffsets[i + 1] - modelInfo.seqOffsets[i];
				alphaAndBetaPass(i, len);

				writeMarginals(i, len);

				scoreAlphaBeta(i, len);

				Arrays.fill(scoreFeatureProductExpectation, 0.0);
				double thisResult = scoreFeatureExpections(i, len, scoreFeatureProductExpectation);
				// Combine the various terms to update the gradient.
				for(int j = 0; j<modelInfo.nFeatures; ++j) {
					grad[j] += scoreFeatureProductExpectation[j] - thisResult * expects[j]; 
					Assert.a(!Double.isNaN(grad[j]));
				}
				if(debug) {
					log.debug(String.format("Iter: %d Seq: %d Expected Score: %g Grad: %s Expected Features: %s Expected Product: %s", iter, i, 
							thisResult, ColtUtil.format(grad), ColtUtil.format(expects), ColtUtil.format(scoreFeatureProductExpectation)));
				}
				result += thisResult;
			}

			// Normalize by the number of positions
			for(int j = 0; j<modelInfo.nFeatures; ++j) {
				grad[j] = grad[j]/modelInfo.totalPositions; 
			}
			result = result/modelInfo.totalPositions; 
			Assert.a(!Double.isNaN(result));
			
			if(debug) {
				log.debug(String.format("Iter: %d Val: %g Grad: %s Weights: %s", iter, result, ColtUtil.format(grad), ColtUtil.format(lambda)));
			}
			else if(log.isInfoEnabled()) {
				log.debug(String.format("Iter: %d Val: %g Grad: %s", iter, result, ColtUtil.format(grad)));
			}

			iter += 1;
		}
		finally {
			logs.close();
		}
		return result;
	}

	private final double scoreFeatureExpections(int seqNum, int len, double[] scoreFeatureProductExpectation) {
		// Now we need to use those alpha and beta vectors to compute the expectations for product of score and features
		// In the process we also need to compute the result
		double result = 0.0;
		int seqOffset = modelInfo.seqOffsets[seqNum];

		scorePosZero(seqNum, scoreFeatureProductExpectation);

		// Put an empty entry in the lookback so the first base has 0's initialized.
		double[] stableState = nextBuffer.stableState;
		Arrays.fill(stableState, 0.0);
		nextBuffer = lookbackBuffer.addFirst(nextBuffer);
		
		int prevPos = 0;
		for(int pos = 1; pos < len; ++pos) {
			int overallPosition = seqOffset + pos;
			cacheProcessor.evaluatePosition(seqNum, pos);

			cacheMi(seqNum, nextBuffer.mi, stableState, nextBuffer.stableState, pos);
			stableState = nextBuffer.stableState;
			nextBuffer = lookbackBuffer.addFirst(nextBuffer);
		
			// Handle the markov features & update the expected value of the similarity score
			FeatureEvaluation nodeEvals = null;
			for(short potential : modelInfo.orderedPotentials) {
				if(potential < modelInfo.nStates) {
					nodeEvals = evals[potential];
				}
				else {
					// This is an edge potential
					int trans = potential - modelInfo.nStates;
					int yprev = modelInfo.transitionFrom[trans];
					int y = modelInfo.transitionTo[trans];
					double ep = allEdgeProb[trans][prevPos]; 
					double npPrev = allNodeProb[yprev][prevPos]; 
					double np = allNodeProb[y][pos]; 
					if(np > 0 && npPrev > 0 && ep >= 0) {
						double scoreExpect = ep * localScoreTransitionCache[trans][overallPosition];	
						result += scoreExpect;

						// Computing marginals from node and edge probabilities is not really safe.  Clamp to a valid probability.
						double conditionalIn = Math.min(1.0, ep / npPrev);
						double conditionalOut = Math.min(1.0, ep / np);

						double inner = scoreExpect;
						inner += conditionalIn*scoreAlpha[yprev][prevPos];
						inner += conditionalOut*scoreBeta[y][pos];

						if(inner == 0.0)
							continue;
						
						// Edges
						FeatureEvaluation potEvals = evals[potential];
						short[] indices = potEvals.index;
						float[] vals = potEvals.value;

						int fCount = 0;
						short index = indices[fCount];
						while(index != -1) {
							// FeatureValue * Score * Marginal Prob
							scoreFeatureProductExpectation[index] += inner*vals[fCount]; 
							if(logs.expectedProductWriter != null) {
								FileUtil.safeWrite(logs.expectedProductWriter, String.format("Seq: %d Pos: %d Edge: %d-%d\tFeat: %d = %g = %g + Val: %g * (s: %g * ep: %g + a: %g * ms: %g + b: %g * me: %g)\n",
										seqNum, pos, yprev, y, index, scoreFeatureProductExpectation[index], scoreFeatureProductExpectation[index]-inner*vals[fCount], vals[fCount],
										localScoreTransitionCache[trans][overallPosition], ep, scoreAlpha[yprev][prevPos], conditionalIn, 
										scoreBeta[y][pos], conditionalOut));
							}
							index = indices[++fCount];
						}

						// Nodes
						indices = nodeEvals.index;
						vals = nodeEvals.value;

						fCount = 0;
						index = indices[fCount];
						while(index != -1) {
							// FeatureValue * Score * Marginal Prob
							scoreFeatureProductExpectation[index] += inner*vals[fCount]; 
							if(logs.expectedProductWriter != null) {
								FileUtil.safeWrite(logs.expectedProductWriter, String.format("Seq: %d Pos: %d State: %d-%d\tNode Feat: %d = %g = %g + Val: %g * (s: %g * ep: %g + a: %g * ms: %g + b: %g * me: %g)\n",
										seqNum, pos, yprev, y, index, scoreFeatureProductExpectation[index], scoreFeatureProductExpectation[index]-inner*vals[fCount], vals[fCount], 
										localScoreTransitionCache[trans][overallPosition], ep, scoreAlpha[yprev][prevPos], conditionalIn, 
										scoreBeta[y][pos], conditionalOut));
							}
							index = indices[++fCount];
						}
					}
				}					
			}

			// Handle the semi-markov features - We need to calculate the probability of each segment
			if(semiMarkov) {
				cacheProcessor.evaluateSegmentsEndingAt(seqNum, pos);
	
				int nSemiMarkovStates = modelInfo.statesWithLookback.length;
				for(int i=0; i<nSemiMarkovStates; ++i) {
					LengthFeatureEvaluation[] lookbacksForState = lengthEvals[i];
					CacheProcessor.StatePotentials statePotentials = modelInfo.statesWithLookback[i];
					byte toNode = statePotentials.state;
					
					int lbIndex=0;
					LengthFeatureEvaluation lengthEval = lookbacksForState[lbIndex];
					int lookback = lengthEval.lookback;
					while(lookback != -1) {
						int beginPos = pos - lookback - 1;
						Assert.a(lengthEval.edgeEvals == null);
	
						double stableScore = localScoreStableCache[toNode][overallPosition] - localScoreStableCache[toNode][overallPosition-lookback];
						double beta = betas[toNode][pos];
						int betaNorm = betaNorms[pos];
	
						// For speed I hand inline RecyclingBuffer.get
						LookbackBuffer segBegin = lookbackBuffer.array[(lookbackBuffer.currentStart+lookback)%lookbackBuffer.length];
						double stableValue = stableState[toNode] - segBegin.stableState[toNode];
	
						// Add in the length based features
						// Handle evaluation of the node potentials
						double nodeValue = stableValue;
						nodeEvals = lengthEval.nodeEval;
						short[] indices = nodeEvals.index;
						float[] vals = nodeEvals.value;
						int ix = 0;
						short index = indices[ix];
						while(index >= 0) {
							nodeValue += vals[ix] * lambda[index];
							index = indices[++ix];
						}
						
						// Look at all the transitions, calculate an inner value for each and sum.  We multiply the whole some
						// by our observed node length features
						double inner = 0.0;
						if(beginPos == -1) {
							double segProb = beta*zInv*exp(starterAlpha[toNode] + nodeValue + NORM_FACTOR*(betaNorm - zNorm));
							if(betaNorm == Integer.MIN_VALUE)
								segProb = 0.0;
	
							inner += segProb * stableScore;
							double segEndMarg = allNodeProb[toNode][pos]-allEdgeProb[modelInfo.selfTransitions[toNode]][pos];
							double conditionalOut = (segEndMarg > 0) ? minMax(segProb/segEndMarg) : 0.0;
							inner += conditionalOut*semiMarkovScoreBeta[toNode][pos];
							if(logs.expectedProductWriter != null) {
								FileUtil.safeWrite(logs.expectedProductWriter, String.format("Seq: %d Pos: 0-%d State: %d\t Inner %g = cOut: %g * sb: %g + stabScore: %g * p: %g "+
										"(b: %g * zInv %g * exp(nodeL:%g + stable: %g + FACTOR*(normB: %d + normZ: %d)))\n",
										seqNum, pos, toNode, inner, conditionalOut, semiMarkovScoreBeta[toNode][pos], stableScore, segProb,
										beta, zInv, nodeValue-stableValue, stableValue+starterAlpha[toNode], betaNorm, zNorm));
							}
						}
						else {
							int nEdges = statePotentials.potentials.length;
							for(int edgeIx=0; edgeIx < nEdges; ++edgeIx) {
								int potential = statePotentials.potentials[edgeIx];
								int trans = potential - modelInfo.nStates;
								int fromNode = modelInfo.transitionFrom[trans];
	
								Assert.a(lengthEval.edgeEvals == null, "Explicit length edge features not supported.");
								
								// Skip semi-Markov self transitions
								if(fromNode == toNode)
									continue;
								
								int invalidIndex = (seqOffset+beginPos+1)*modelInfo.nPotentials;
								if(invalidTransitions[invalidIndex + potential]) {
									continue;
								}
	
								double prevAlpha = alphas[beginPos][fromNode];
								int prevAlphaNorm = alphaNorms[beginPos];
								
								double transitionValue = segBegin.mi[trans];
								double segProb = prevAlpha * beta * zInv * exp(nodeValue + transitionValue + NORM_FACTOR*(prevAlphaNorm + betaNorm - zNorm));
								if(prevAlphaNorm == Integer.MIN_VALUE || betaNorm == Integer.MIN_VALUE)
									segProb = 0.0;
								if(Double.isNaN(segProb))
									Assert.a(false, "Bad Segment Prob. Seq ", seqNum, " Pos ",prevPos, "-", pos);
	
								double segmentScore = localScoreTransitionCache[trans][seqOffset + beginPos+1] + stableScore;
								inner += segProb * segmentScore;

								double prevSegMarg = allNodeProb[fromNode][beginPos] - ((modelInfo.maxStateLengths[fromNode]>1) ?
										allEdgeProb[modelInfo.selfTransitions[fromNode]][beginPos] : 0);
								double conditionalIn = (prevSegMarg > 0) ? minMax(segProb/prevSegMarg) : 0.0;
								inner += conditionalIn*semiMarkovScoreAlpha[fromNode][beginPos];

								double segEndMarg = allNodeProb[toNode][pos]-allEdgeProb[modelInfo.selfTransitions[toNode]][pos];
								double conditionalOut = (segEndMarg > 0) ? minMax(segProb/segEndMarg) : 0.0;
								inner += conditionalOut*semiMarkovScoreBeta[toNode][pos];
								if(logs.expectedProductWriter != null) {
									FileUtil.safeWrite(logs.expectedProductWriter, String.format("Seq: %d Pos: %d-%d State: %d-%d\t Inner %g = cIn: %g * sa: %g + cOut: %g * sb: %g + score: %g (trans: %g + stab: %s) * p: %g "+
											"(a: %g * b: %g * zInv %g * exp(nodeL:%g + stable: %g + trans: %g + FACTOR*(normA: %d + normB: %d + normZ: %d)))\n",
											seqNum, beginPos+1, pos, fromNode, toNode, inner, conditionalIn, semiMarkovScoreAlpha[fromNode][beginPos], conditionalOut, semiMarkovScoreBeta[toNode][pos], segmentScore, localScoreTransitionCache[trans][seqOffset + beginPos+1], stableScore, segProb,
											prevAlpha, beta, zInv, nodeValue-stableValue, stableValue, transitionValue, prevAlphaNorm, betaNorm, zNorm));
								}
							}
						}
						                                           					
						// Once the inner value has been computed, multiply it by all of the observed feature values.
						FeatureEvaluation lengthNodeEvals = lengthEval.nodeEval;
						indices = lengthNodeEvals.index;
						vals = lengthNodeEvals.value;
						
						ix = 0;
						index = indices[ix];
						while(index >= 0) {
							scoreFeatureProductExpectation[index] += inner*vals[ix]; 
							if(logs.expectedProductWriter != null) {
								FileUtil.safeWrite(logs.expectedProductWriter, String.format("Seq: %d Pos: %d-%d State: %d\tFeat: %d = %g = %g + Val: %g * Inner: %g\n",
										seqNum, beginPos+1, pos, toNode, index, scoreFeatureProductExpectation[index], scoreFeatureProductExpectation[index]-inner*vals[ix], vals[ix], inner));
							}
							index = indices[++ix];
						}
						Assert.a(lengthEval.edgeEvals == null, "Explicit length edges are not supported.");
						
						++lbIndex;
						lengthEval = lookbacksForState[lbIndex];
						lookback = lengthEval.lookback;
					}
				}
			}
			prevPos = pos;
		}
		return result;
	}		
		
	void scorePosZero(int seqNum, double[] scoreFeatureProductExpectation) {
		cacheProcessor.evaluatePosition(seqNum, 0);
		for(int state = 0; state < modelInfo.nStates; ++state) {
			// This is a node potential
			double inner = scoreBeta[state][0];
			FeatureEvaluation potEvals = evals[state];
			short[] indices = potEvals.index;
			float[] vals = potEvals.value;

			int fCount = 0;
			short index = indices[fCount];
			while(index != -1) {
				// FeatureValue * Score * Marginal Prob
				scoreFeatureProductExpectation[index] += inner*vals[fCount]; 
				if(logs.expectedProductWriter != null) {
					FileUtil.safeWrite(logs.expectedProductWriter, String.format("Seq: %d Pos: 0 State: %d\tFeat: %d = %g = %g + Val: %g * Beta[%d][%d]: %g:\n",
							seqNum, state, index, scoreFeatureProductExpectation[index], scoreFeatureProductExpectation[index]-inner*vals[fCount], vals[fCount], 0, state, inner));
				}
				index = indices[++fCount];
			}
		}

		// Handle the semi-markov features - We avoid having to calculate segment probabilities
		if(semiMarkov) {
			cacheProcessor.evaluateSegmentsEndingAt(seqNum, 0);
	
			int nSemiMarkovStates = modelInfo.statesWithLookback.length;
			for(int i=0; i<nSemiMarkovStates; ++i) {
				LengthFeatureEvaluation[] lookbacksForState = lengthEvals[i];
				LengthFeatureEvaluation lengthEval = lookbacksForState[0];
				if(lengthEval.lookback != -1) {
					Assert.a(lengthEval.lookback == 0);
					CacheProcessor.StatePotentials statePotentials = modelInfo.statesWithLookback[i];
					byte state = statePotentials.state;
	
					double inner = semiMarkovScoreBeta[state][0];
					short[] indices = lengthEval.nodeEval.index;
					float[] vals = lengthEval.nodeEval.value;
					int fCount = 0;
					short index = indices[fCount];
					while(index != -1) {
						// FeatureValue * Score * Marginal Prob
						scoreFeatureProductExpectation[index] += inner*vals[fCount]; 
						if(logs.expectedProductWriter != null) {
							FileUtil.safeWrite(logs.expectedProductWriter, String.format("Seq: %d Pos: 0 State: %d\tLen. Feat: %d = %g = %g + Val: %g * Beta[%d][%d]: %g:\n",
									seqNum, state, index, scoreFeatureProductExpectation[index], scoreFeatureProductExpectation[index]-inner*vals[fCount], vals[fCount], 0, state, inner));
						}
						index = indices[++fCount];
					}
					Assert.a(lookbacksForState[1].lookback == -1);
				}
			}
		}
	}

	private final void scoreAlphaBeta(int seqNum, int len) {
		int seqOffset = modelInfo.seqOffsets[seqNum];
		
		// Do another backward & forward pass to compute the score alpha & betas
		// scoreBeta is being defined from 0 to len-1
		int pos = len-1;
		int prevPos;
		for (int y=0; y<modelInfo.nStates; y++) {
			Arrays.fill(scoreBeta[y], 0.0);
			if(semiMarkov)
				Arrays.fill(semiMarkovScoreBeta[y], 0.0);
		}
		for (prevPos = len-2; prevPos >= 0; --prevPos) {
			for (int trans=0; trans<modelInfo.nTransitions; trans++) {
				int yprev = modelInfo.transitionFrom[trans];
				int y = modelInfo.transitionTo[trans];
				double ep = allEdgeProb[trans][prevPos];
				double np = allNodeProb[y][pos];
				if(np > 0 && ep >= 0) {
					// Computing marginals from node and edge probabilities is not really safe.  Clamp to a valid probability.
					double conditional = Math.min(1.0, ep / np);
					double update = ep*localScoreTransitionCache[trans][seqOffset+pos] + conditional*scoreBeta[y][pos];
					scoreBeta[yprev][prevPos] += update;
					if(modelInfo.maxStateLengths[y]>1 && y != yprev)
						semiMarkovScoreBeta[yprev][prevPos] += update;
				}					
			}
			pos = prevPos;
		}
		
		// scoreAlpha is being defined from 0 to len-2 (we never use the last alpha position)
		prevPos = 0;
		for (int y=0; y<modelInfo.nStates; y++) {
			Arrays.fill(scoreAlpha[y], 0.0);
			if(semiMarkov)
				Arrays.fill(semiMarkovScoreAlpha[y], 0.0);
		}
		for (pos=1; pos<len-1; pos++) {
			for (int trans=0; trans<modelInfo.nTransitions; trans++) {
				int yprev = modelInfo.transitionFrom[trans];
				int y = modelInfo.transitionTo[trans];
				double ep = allEdgeProb[trans][prevPos];
				double np = allNodeProb[yprev][prevPos];
				if(np > 0 && ep >= 0) {
					// Computing marginals from node and edge probabilities is not really safe.  Clamp to a valid probability.
					double conditional = Math.min(1.0, ep / np);

					double update = conditional*scoreAlpha[yprev][prevPos];
					update += ep*localScoreTransitionCache[trans][seqOffset+pos];
					scoreAlpha[y][pos] += update;
					if(modelInfo.maxStateLengths[y]>1) {
						double nodeMarg = allNodeProb[y][pos];
						if(nodeMarg > 0) {
							double outConditional = minMax((nodeMarg-allEdgeProb[modelInfo.selfTransitions[y]][pos])/nodeMarg);
							semiMarkovScoreAlpha[y][pos] += update*outConditional;
						}
					}
						
					if(logs.scoreAlphaWriter != null) {
						FileUtil.safeWrite(logs.scoreAlphaWriter, String.format("Seq: %d alpha[%d][%d] = %g = %g + Pr: %g * alpha[%d][%d] %g + Pr: %g * Score: %g\n",
								seqNum, pos, y, scoreAlpha[y][pos], scoreAlpha[y][pos]-update, ep/allNodeProb[yprev][prevPos], prevPos, yprev, scoreAlpha[yprev][prevPos], ep, localScoreTransitionCache[trans][seqOffset+pos]));
					}
				}
			}	
			prevPos = pos;
		}
	}
	
	void fillScoreCache(List<? extends TrainingSequence<?>> data) {
		localScoreStableCache = new double[modelInfo.nStates][modelInfo.totalPositions];
		localScoreTransitionCache = new double[modelInfo.nTransitions][modelInfo.totalPositions];
		
		int overallPosition = 0;
		for(int i=0; i<data.size(); ++i) {
			TrainingSequence seq = data.get(i);
			++overallPosition;
			for(int pos= 1; pos<seq.length(); ++pos ) {
				for(int transition = 0; transition < modelInfo.nTransitions; ++transition) {
					int from = modelInfo.transitionFrom[transition];
					int to = modelInfo.transitionTo[transition];
					double localScore = score.evaluate(from, to, seq, pos);
					localScoreTransitionCache[transition][overallPosition] = localScore;
					if(from == to)
						localScoreStableCache[to][overallPosition] = localScoreStableCache[to][overallPosition-1] + localScore;
				}
				++overallPosition;
			}
		}
	}
	
	private static final double minMax(final double val) {
		double ret = (val <= 0.0) ? 0.0 : ((val >= 1.0) ? 1.0 : val);
		if(Double.isNaN(ret)) {
			Assert.a(false, "Min max called on "+val);
		}
		return ret;
	}
	
	private final void writeMarginals(int i, int len) {
		if(logs.marginalsWriter != null) {
			for(int pos=0; pos < len; ++pos) {
				FileUtil.safeWrite(logs.marginalsWriter, String.format("Seq %d Pos %d -", i, pos));
				for(short potential : modelInfo.orderedPotentials) {
					if(potential < modelInfo.nStates) {
						FileUtil.safeWrite(logs.marginalsWriter, String.format(" State %d: %e", potential, allNodeProb[potential][pos]));
					}
					else {
						int trans = potential - modelInfo.nStates;
						FileUtil.safeWrite(logs.marginalsWriter, String.format(" Edge %d-%d: %e", modelInfo.transitionFrom[trans], modelInfo.transitionTo[trans], allEdgeProb[trans][pos]));
					}
				}
				FileUtil.safeWrite(logs.marginalsWriter, "\n");
			}
		}
	}
	
	/** gets the local similarity score function used to score each position in every path.
	 * @return the configured score function
	 */
	public LocalPathSimilarityScore getScore() {
		return score;
	}

	/** sets the local similarity score function used to score each position in every path.  This
	 * is usually specified in the XML configuration file.
	 * @param score the score function to use
	 */
	public void setScore(LocalPathSimilarityScore score) {
		this.score = score;
	}

	public String getScoreAlphaFile() {
		return logs.scoreAlphaFile;
	}

	public void setScoreAlphaFile(String scoreAlphaFile) {
		logs.scoreAlphaFile = scoreAlphaFile;
	}

	public String getExpectedProductFile() {
		return logs.expectedProductFile;
	}

	public void setExpectedProductFile(String expectedProductFile) {
		logs.expectedProductFile = expectedProductFile;
	}

	public String getMarginalsFile() {
		return logs.marginalsFile;
	}
	
	public void setMarginalsFile(String marginalsFile) {
		logs.marginalsFile = marginalsFile;
	}
}
