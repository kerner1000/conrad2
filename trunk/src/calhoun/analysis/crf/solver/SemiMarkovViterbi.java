package calhoun.analysis.crf.solver;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.CRFInference;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.analysis.crf.solver.CacheProcessor.FeatureEvaluation;
import calhoun.analysis.crf.solver.CacheProcessor.LengthFeatureEvaluation;
import calhoun.analysis.crf.solver.CacheProcessor.SolverSetup;
import calhoun.util.Assert;
import calhoun.util.ColtUtil;

/** viterbi algorithm for semi-Markov CRFs.  Uses a {@link CacheProcessor} to handle the complexities of evaluation.
 * This is similar to a markov viterbi search, except we have to search over all valid segments to fill in one entry 
 * in the matrix instead of just looking at the last position.
 */
public class SemiMarkovViterbi implements CRFInference {
	private static final Log log = LogFactory.getLog(SemiMarkovViterbi.class);
	boolean debug = log.isDebugEnabled();

	double[] lambda;

	private double[] bestScore;
	private int[] backPointers;
	int nStates;
	
	SolverSetup modelInfo;
	CacheProcessor cacheProcessor;
	FeatureEvaluation[] evals;
	LengthFeatureEvaluation[][] lengthEvals;
	boolean[] invalidTransitions;
	int[] selfTransitions;

	RecyclingBuffer<double[]> stableStates;
	double[] stableVector;
	
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

	public InferenceResult predict(ModelManager fm, InputSequence<?> seq, double[] lambda) {
		this.lambda = lambda;

		cacheProcessor.setInputData(fm, seq);
		modelInfo = cacheProcessor.getSolverSetup();
		nStates = modelInfo.nStates;
		Assert.a(modelInfo.maxStateLengths.length == nStates, "Maximum state lengths array was length ("+modelInfo.maxStateLengths.length+").  Must have one entry for each state "+modelInfo.nStates+")");
		evals = cacheProcessor.getFeatureEvaluations();
		lengthEvals = cacheProcessor.getLengthFeatureEvaluations();
		invalidTransitions = cacheProcessor.getInvalidTransitions();
		
		int len = seq.length();
		
		selfTransitions = new int[nStates];
		for(int i=0; i<nStates; ++i) {
			selfTransitions[i] = modelInfo.selfTransitions[i];
		}

		// Circular buffers of the previous mi matrices.
		RecyclingBuffer<double[]> mis = new RecyclingBuffer<double[]>(new double[modelInfo.maxLookback][modelInfo.nTransitions]);
		double[] nextMi = new double[modelInfo.nTransitions];

		// Circular buffers of the values of staying in a stable state for a given period of time
		stableStates = new RecyclingBuffer<double[]>(new double[modelInfo.maxLookback][nStates]);
		stableVector = new double[nStates];
		
		bestScore = new double[len*nStates];
		backPointers = new int[len*nStates];
		int[] backLengths = new int[len*nStates];
		
		for (int pos = 0; pos < len; pos++) {
			/* compute weighted features.  These are for transitions at the current base for non-length dependent features. */
						
			// Save the first Ri matrix, since it has the initial probabilies
			if(pos == 0) {
				computeSparseMi(seq, pos, null, stableVector);
				stableVector = stableStates.addFirst(stableVector);
			}
			else {
				computeSparseMi(seq, pos, nextMi, null);

				// Use this transition matrix to update the stable vectors 
				updateStableBuffer(nextMi);				

				// Add this into the list of saved Mi matrices
				nextMi = mis.addFirst(nextMi);				
			}				

			double[] latestStable = stableStates.get(0);
			double[] latestMi = mis.get(0);

			// Do states without lookback first.
			for(CacheProcessor.StatePotentials potentials : modelInfo.statesWithoutLookback) {
				int state = potentials.state;
				double max = Double.NEGATIVE_INFINITY;
				int invalidIndex = pos*modelInfo.nPotentials;
				int bestLookback = 0;
				int bestPrevState = -2;
				
				// Check that the node is valid, otherwise put in a dummy entry. 
				if(!invalidTransitions[invalidIndex + state]) {
					if(pos == 0) {
						if(debug) 
							log.debug(String.format("Pos: %d State: %d %.2f", pos, state, latestStable[state]));
						// For the first position, we just store the starting potential
						max = latestStable[state];
						bestPrevState = -1;
					}
					else {
						for(byte edgePotential : potentials.potentials) {
							if(invalidTransitions[invalidIndex + state])
								continue;

							int transition = edgePotential - nStates;
							int prevState = modelInfo.transitionFrom[transition];
							double transitionCost = latestMi[transition];

							if(Double.isInfinite(transitionCost)) {
								// This transition was invalid at this location, ignore.
								continue;
							}

							double previous = bestScore[nStates *(pos-1) + prevState]; 
							double current = previous + transitionCost;
							if(debug)
								log.debug(String.format("Pos: %d Trans: %d-%d %.2f (Prev: %.2f + Trans: %.2f)", pos, prevState, state, current, previous, transitionCost));
							if(current > max) {
								max = current;
								bestPrevState = prevState;
							}
						}
					}
				}
				
				// Fill in our best entry
				int index= pos*nStates + state;
				//log.debug(String.format("bestScore[%d] = %.2f", index, max));
				bestScore[index] = max; 
				backPointers[index] = bestPrevState;
				backLengths[index] = bestLookback+1;
			}
			
			// Now repeat for states with lookback.
			for(int i=0; i<modelInfo.statesWithLookback.length; ++i) {
				CacheProcessor.StatePotentials potentials = modelInfo.statesWithLookback[i];
				LengthFeatureEvaluation[] lookbacksForState = lengthEvals[i];

				int state = potentials.state;
				double max = Double.NEGATIVE_INFINITY;
				int bestLookback = -1;
				int bestPrevState = -2;
				
				cacheProcessor.evaluateSegmentsEndingAt(0, pos);

				int lbIndex=0;
				LengthFeatureEvaluation lengthEval = lookbacksForState[lbIndex];
				int lookback = lengthEval.lookback;
				while(lookback != -1) {
					//log.info("Pos: "+pos+"\t State: "+modelInfo.statesWithLookback[i].state+"\t Lookback: "+lookback);

					double[] lookbackMi = mis.get(lookback); 
					double[] lookbackStable = stableStates.get(lookback);

					// Handle evaluation of the node potentials
					FeatureEvaluation nodeEvals = lengthEval.nodeEval;
					short[] indices = nodeEvals.index;
					float[] vals = nodeEvals.value;
					int ix = 0;
					short index = indices[ix];
					double nodePotential = 0.0;
					while(index >= 0) {
						nodePotential += vals[ix] * lambda[index];
						index = indices[++ix];
					}
					Assert.a(index != Short.MIN_VALUE, "Node lengths should only be returned in the cache if they are valid");

					int prevPos = pos - lookback - 1;
					if(prevPos < 0) {
						// Examine the case where this is the first segment in the sequence.
						// Score is the sum of non-length dependent fetaures plus length features.
						double current = latestStable[state] + nodePotential;
						if(debug)
							log.debug(String.format("Pos: %d Lb: %d State: %d %.2f (Stable: %.2f + Node: %.2f)", pos, lookback, state, current, latestStable[state], nodePotential));
						
						// Check if this is our best so far
						if(current > max) {
							max = current;
							bestLookback = lookback;
							bestPrevState = -1;
						}
					}
					else {
						// If this is not the first segment, we need to deal with edges coming into this segment
						FeatureEvaluation[] edgeEvals = lengthEval.edgeEvals;
						int nEdges = potentials.potentials.length;
						for(int edgeIx=0; edgeIx < nEdges; ++edgeIx) {
							int potential = potentials.potentials[edgeIx];
							int trans = potential - modelInfo.nStates;
							int fromNode = modelInfo.transitionFrom[trans];
							// Skip semi-Markov self transitions
							if(fromNode == state)
								continue;

							double edgeVal = 0.0;

							if(edgeEvals == null) {
								// If the cache processor does not have edge evaluations
								// Just check if this transition is legal based on the invalid transitions matrix
								if(invalidTransitions[(prevPos+1)*modelInfo.nPotentials + potential]) {
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
								if(index == Short.MIN_VALUE) {
									log.info("SHORT.MIN_VALUE");
									continue;
								}
								while(index != -1) {
									edgeVal += vals[ix] * lambda[index];
									index = indices[++ix];
								}
							}
							
							double prevBest = bestScore[nStates *(pos-(lookback+1)) + fromNode];
							double stable = latestStable[state]-lookbackStable[state];
							
							// Renormalize and update the exp value.
							double current = prevBest + nodePotential + edgeVal + stable + lookbackMi[trans];
							if(debug)
								log.debug(String.format("Pos: %d Lb: %d Trans: %d-%d %.4f (Prev: %.4f + Stable: %.4f + Trans: %.4f + Node: %.4f + Edge: %.4f)", pos, lookback, fromNode, state, 
										current, prevBest, stable, lookbackMi[trans], nodePotential, edgeVal));
							
							if(current == Double.NEGATIVE_INFINITY) {
								// The state was not valid at the previous position
								continue;
							}

							if(current > max) {
								//log.info("Selected Pos: "+pos+" Edge: "+prevState+"-"+state+" "+current+" vs. "+max+" Prev. "+bestScore[nStates *(pos-(lookback+1)) + prevState]+" Trans. "+transitionCost);
								max = current;
								bestLookback = lookback;
								bestPrevState = fromNode;
							}
						}
					}
					++lbIndex;
					lengthEval = lookbacksForState[lbIndex];
					lookback = lengthEval.lookback;
				}
				
				// Fill in our best entry
				//log.info(String.format("Pos: %d State: %d BestScore: "+max+" BackPointer: %d", pos, state, bestPrevState));
				int index= pos*nStates + state; 
				bestScore[index] = max; 
				backPointers[index] = bestPrevState;
				backLengths[index] = bestLookback+1;
			}
		}
		//printForwardPass(bestScore, backPointers, backLengths, nStates, seq.length());
		//log.info(ColtUtil.format(bestScore));
		//log.info(StringUtils.join(backPointers," "));
		//log.info(backLengths);
				
		// Now that we have the matrix, trace back to get the best path.
		int[] ret = new int[len];
		int pos = len-1;
		int state = ColtUtil.maxInColumn(bestScore, nStates, len-1);
		Assert.a(state != -2, "No valid paths");

		while(pos >= 0) {
			int stateLen = backLengths[pos*nStates + state];
			int prevState = backPointers[pos*nStates + state];
			//log.info(String.format("State: %d, Len: %d, Ends At: %d", state, stateLen, pos));
			for(int i = 0; i < stateLen; ++i) {
				ret[pos] = state;
				pos--;
			}
			state = prevState; 
		}		
		Assert.a(pos == -1);
		InferenceResult inferenceResult = new InferenceResult();
		inferenceResult.hiddenStates = ret;
		inferenceResult.bestScores = new double[nStates];
		System.arraycopy(bestScore, nStates*(len-1), inferenceResult.bestScores, 0, nStates );
		return inferenceResult;
	}

	/** Updates the stableStates list based on the current position.  The stableStates list contains the 
	 cost of the non-length dependent features for the duration of the lookback. 
	 */
	void updateStableBuffer(double[] nextMi) {
		// Recycle matrices at the end
		double[] prevState = stableStates.get(0);
		for(int ix = 0; ix < nStates; ++ix) {
			if(modelInfo.maxStateLengths[ix] > 1) {
				int trans = selfTransitions[ix];
				if(trans != -1) {
					if (Double.isInfinite(nextMi[trans]))
						stableVector[ix] = prevState[ix];
					else
						stableVector[ix] = prevState[ix] + nextMi[trans];				
					//log.info(String.format("stableState[%d] = %f = %f + %f", ix, stableState[ix], prevState[ix], nextMi[trans]));
				}
			}
		}

		// Now add this to the beginning.
		stableVector = stableStates.addFirst(stableVector);
	}

	void computeSparseMi(InputSequence seq, int pos, double[] mi, double[] ri) {
		cacheProcessor.evaluatePosition(0, pos);
		double nodeVal = Double.NaN;
		int invalidIndex = pos*modelInfo.nPotentials;
		for(short potential : modelInfo.orderedPotentials) {
			boolean invalid = invalidTransitions[invalidIndex + potential];
			double features = invalid ? Double.NEGATIVE_INFINITY : 0.0;

			// Add up all features for this potential.
			FeatureEvaluation potEvals = evals[potential];
			short[] indices = potEvals.index;
			float[] vals = potEvals.value;
			int i = 0;
			short index = indices[i];
			while(index != -1) {
				// An invalid potential is indicated by a feature value of Short.MAX_VALUE
				features += (index == Short.MIN_VALUE) ? Double.NEGATIVE_INFINITY : vals[i]*lambda[index]; 
				index = indices[++i];
			}
			
			if(potential < modelInfo.nStates) {
				nodeVal = features;
				if(ri != null) {
					ri[potential] = nodeVal;
				}
			}
			else {
				//log.debug(String.format("Mi[%d, %d] = %f, adding in %f to get %f", feat.yprev(), feat.y(), val, feat.value(), val*exp(feat.value()*param[feat.index()])));
				int transition = potential - modelInfo.nStates;
				double val = features + nodeVal;
				if(mi != null)
					mi[transition] = val;
			}
		}		
	}
}
