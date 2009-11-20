package calhoun.analysis.crf.solver;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.CRFInference;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.SemiMarkovSetup;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.analysis.crf.solver.check.FeatureCalculator;
import calhoun.analysis.crf.solver.check.TransitionInfo;
import calhoun.util.Assert;
import calhoun.util.ColtUtil;
import calhoun.util.DenseIntMatrix2D;

/** viterbi algorithm for semi-Markov CRFs.  Does not perform any caching, and <b>does not correctly handle
 * NodeBoundary features</b>.  Use only if you have no NodeBoundary features and are running into memory limitations.
 * Requires that a {@link SemiMarkovSetup} be configured.
 */
public class SemiMarkovViterbiNoCache implements CRFInference {
	private static final Log log = LogFactory.getLog(SemiMarkovViterbi.class);
	boolean debug = log.isDebugEnabled();

	private boolean allPaths;
	private double[] bestScore;
	private int[] backPointers;
	DenseIntMatrix2D backLengths;
	short[] maxStateLengths;
	boolean ignoreSemiMarkovSelfTransitions;
	
	int maxLookback = 1;
	int nStates;

	TransitionInfo transitions;
	int[] selfTransitions;
	
	class LengthTransitionInfo extends TransitionInfo {
		short[] localMaxStateLengths;
		LengthTransitionInfo(ModelManager fm, boolean allPaths, short[] maxStateLengths, boolean ignoreSemiMarkovSelfTransitions) {
			this.localMaxStateLengths = maxStateLengths;
			this.ignoreSemiMarkovSelf = ignoreSemiMarkovSelfTransitions;
			initTrans(fm, allPaths);
		}
		@Override
		protected boolean allowSelf(int state) {
			return localMaxStateLengths[state] > 1;
		}
	}
	
	/** true if all paths (valid and invalid) are to be evaluated during the viterbi search.  Defaults to false.
	 * @return true if all paths are to be examined
	 */
	public boolean isAllPaths() {
		return allPaths;
	}

	/** sets whether all paths (valid and invalid) are to be evaluated during the viterbi search.  Defaults to false.
	 * @param allPaths allPath true if all paths are to be examined
	 */
	public void setAllPaths(boolean allPaths) {
		this.allPaths = allPaths;
	}

	/** sets up the semiMarkov parameters, such as maximum segment lengths.  These should be the same as the model was trained with.
	 * @param setup the parameters to use during the semi-Markov viterbi search
	 */
	public void setSemiMarkovSetup(SemiMarkovSetup setup) {
		maxStateLengths = setup.getMaxLengths();
		ignoreSemiMarkovSelfTransitions = setup.isIgnoreSemiMarkovSelfTransitions();
	}

	public InferenceResult predict(ModelManager fm, InputSequence<?> seq, double[] lambda) {		
		nStates = fm.getNumStates();
		int len = seq.length();
		if(maxStateLengths == null) {
			log.info("No state lengths set - standard viterbi search.");
			maxStateLengths = new short[nStates];
			Arrays.fill(maxStateLengths, (short) 1);
		}
		
		// Determine the longest lookback we might have to do.
		for(int i : maxStateLengths) {
			maxLookback = Math.max(i, maxLookback);
		}

		transitions = new LengthTransitionInfo(fm, allPaths, maxStateLengths, ignoreSemiMarkovSelfTransitions);
		int[][] statePotentials = getStatePotentials(transitions);
		selfTransitions = new int[nStates];
		for(int i=0; i<nStates; ++i) {
			selfTransitions[i] = transitions.transitionIndex.getQuick(i, i);
		}
		
		FeatureCalculator calc = new FeatureCalculator(fm, lambda, transitions);

		// Buffer of the previous mi matrices.
		LinkedList<double[]> mis = new LinkedList<double[]>();
		// Buffer of the values of staying in a stable state for a given period of time
		LinkedList<double[]> stableStates = new LinkedList<double[]>();
		
		// The next MI matrix we plan to use (may be recycled from the end of the list
		double[] nextMi = new double[transitions.nTransitions];

		double[] Ri = new double[nStates];
		
		bestScore = new double[seq.length()*nStates];
		backPointers = new int[seq.length()*nStates];
		int[] localBackLengths = new int[seq.length()*nStates];
		
		for (int pos = 0; pos < len; pos++) {
			/* compute weighted features.  These are for transitions at the current base for non-length dependent features. */
			calc.computeSparseMi(seq, pos, nextMi, Ri);
						
			// Save the first Ri matrix, since it has the initial probabilies
			if(pos == 0) {
				stableStates.add(Ri);
				Ri = new double[nStates];
			}
			else {
				// Use this transition matrix to update the stable vectors 
				updateStableBuffer(stableStates, nextMi);				

				// Add this into the list of saved Mi matrices
				nextMi = updateMiBuffer(mis, nextMi);				
			}				

			double[] latestStable = stableStates.getFirst();
			// Now investigate each state to determine the optimal path to this point.
			for(int state = 0; state<nStates; ++state) {
				int[] transitionPotentials = statePotentials[state];
				
				// For each state, determine how far we should look back 
				int lookbackSize = maxStateLengths[state];
				
				// Examine each length, starting with 1 and increasing.
				double max = Double.NEGATIVE_INFINITY;
				int bestLookback = -1;
				int bestPrevState = -2;
				
				// As we go through the lookbacks for the previous state, go through the stored Mi matrics and stable positions
				Iterator<double[]> miIter = mis.iterator();
				Iterator<double[]> stableIter = stableStates.iterator();
				for(int lookback = 0; lookback < lookbackSize; ++lookback) {
					// Find the starting position we are evaluating for this state
					int startPos = pos - lookback;
					//Assert.a(startPos >= 0);
					if(startPos == 0) {
						// In this case we have the extra check that we should have the same number of previous Mis as our length
						Assert.a(mis.size() == lookback, "More Mi matrices in history than there are previous positions in the sequence.");
						
						// Examine the case where this is the first position in the sequence.
						// Score is the initial score plus the score of staying in that feature for this long.
						double current = latestStable[state];
						// Now add length dependent features
						calc.result.evaluateNodeLength(seq, pos, lookback+1, state);
						double lengthCost = calc.calcRet(false);
						current += lengthCost;
						
						// Check if this is our best so far
						if(current > max) {
							max = current;
							bestLookback = lookback;
							bestPrevState = -1;
						}
						// Don't bother looking back any farther
						break;
					}
					else {
						double nodeLengthCost = Double.NaN;
						// This state starts after the beginning of the sequences.
						// Need to check all legal transitions from previous states.
						// This includes transitions to self although I'm not sure I like that for states with length models.
						double[] lookbackMi = miIter.next(); 
						double[] lookbackStable = stableIter.next();
						for(int transition : transitionPotentials) {
							int prevState = transitions.transitionFrom[transition];
							if(prevState == state && lookbackSize > 1) {
								// For states with explicit length distributions, ignore transition to self.
								continue;
							}
							double transitionCost = lookbackMi[transition];
							if(Double.isInfinite(transitionCost)) {
								// This transition was invalid at this location, ignore.
								continue;
							}
							if(Double.isNaN(nodeLengthCost)) {
								// Only calculate this if you know you need it.
								nodeLengthCost = calc.calcNodeLengthValue(seq, pos, lookback+1, state);
							}
							double current = bestScore[nStates *(pos-(lookback+1)) + prevState];
							if(current == Double.NEGATIVE_INFINITY) {
								// The state was not valid at the previous position
								continue;
							}
							double stable = latestStable[state]-lookbackStable[state];
							calc.result.evaluateEdgeLength(seq, pos, lookback+1, prevState, state);
							double lengthCost = nodeLengthCost + calc.calcRet(false);
							
							current += transitionCost + stable + lengthCost;
							if(current > max) {
								//log.info("Selected Pos: "+pos+" Edge: "+prevState+"-"+state+" "+current+" vs. "+max+" Prev. "+bestScore[nStates *(pos-(lookback+1)) + prevState]+" Trans. "+transitionCost);
								max = current;
								bestLookback = lookback;
								bestPrevState = prevState;
							}
						}
					}
				}
				
				// Fill in our best entry
				//log.info(String.format("Pos: %d State: %d BestScore: "+max+" BackPointer: %d", pos, state, bestPrevState));
				int index= pos*nStates + state; 
				bestScore[index] = max; 
				backPointers[index] = bestPrevState;
				localBackLengths[index] = bestLookback+1;
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
			int stateLen = localBackLengths[pos*nStates + state];
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
	void updateStableBuffer(LinkedList<double[]> stableStates, double[] nextMi) {
		// Recycle matrices at the end
		double[] stableState;
		if(stableStates.size() > maxLookback) {
			stableState = stableStates.removeLast();
		}
		else {
			stableState = new double[nStates];
		}
		
		double[] prevState = stableStates.getFirst();
		for(int ix = 0; ix < nStates; ++ix) {
			if(maxStateLengths[ix] > 1) {
				int trans = selfTransitions[ix];
				if(trans != -1) {
					if (Double.isInfinite(nextMi[trans]))
						stableState[ix] = prevState[ix];
					else
						stableState[ix] = prevState[ix] + nextMi[trans];				
					//log.info(String.format("stableState[%d] = %f = %f + %f", ix, stableState[ix], prevState[ix], nextMi[trans]));
				}
			}
		}

		// Now add this to the beginning.
		stableStates.addFirst(stableState);
	}

	double[] updateMiBuffer(LinkedList<double[]> mis, double[] nextMi) {
		mis.addFirst(nextMi);
		if(mis.size() > maxLookback) {
			nextMi = mis.removeLast();
		}
		else {
			nextMi = new double[transitions.nTransitions];
		}
		return nextMi;
	}
	
	// Used for debugging only.
	/*
	 * private void printForwardPass(double[] bestScore, int[] backPointers, int[] backLengths, int numStates, int seqLen)
	{
		try {
			//Writer fout = new BufferedWriter(new FileWriter("test/working/crf_forwardPass.txt"));
			Writer fout = new BufferedWriter(new FileWriter("crf_forwardPass.txt"));
			
			int i, pos, st;
			//double[][] bestScores = new double[seqLen][numStates];
			//int[][] bestIndices  = new int[seqLen][numStates];
			
			fout.write("Viterbi Map from Forward Pass\n");
			
			fout.write("\n");

			for (pos=0; pos<seqLen; pos++)
			{
				for (st=0; st<numStates; st++)
				{
					fout.write( (pos*numStates + st) + "\t" + pos + "\t" + st + "\t");
					if (bestScore[pos*numStates + st] == 0)
						fout.write(String.format("%1$11.0f", bestScore[pos*numStates + st]) + "  \t" + "." + "\n");
					else if (bestScore[pos*numStates + st] == Double.NEGATIVE_INFINITY)
						fout.write(String.format("%1$11.0f", 0.0f) + "  \t" +  "." + "\n");
					else {
						//fout.write(String.format("%1$11.2f", bestScore[pos*numStates + st]) + "  \t" + "--" + "\n");

						if (backPointers[pos*numStates + st] == -1)
							fout.write(String.format("%1$11.2f", bestScore[pos*numStates + st]) + "  \t" + "--" + "\t" +  backPointers[pos*numStates + st] + "\n");
						else 
							fout.write(String.format("%1$11.2f", bestScore[pos*numStates + st]) + "  \t" + (pos - backLengths[pos*numStates + st]) + "\t" +  backPointers[pos*numStates + st] + "\n");
				
					}
				}
				//fout.write("\n");
			}
			fout.close();
		
		} catch (IOException e) {
			throw new RuntimeException("Error writing alpha pass");
		}	
	}
	*/

	int[][] getStatePotentials(TransitionInfo transitions) {
		int[][] statePotentials = new int[nStates][];
		int currentState = -1;
		List<Integer> currentList = null;
		for(int potential : transitions.orderedPotentials) {
			if(potential < nStates) {
				if(currentState != -1) {
					statePotentials[currentState] = toIntArray(currentList);
				}
				currentState = potential;
				currentList = new ArrayList<Integer>();
			}
			else {
				currentList.add(potential-nStates);
			}
		}
		statePotentials[currentState] = toIntArray(currentList);
		return statePotentials;
	}
	
	int[] toIntArray(List<Integer> list) {
		int[] ret = new int[list.size()];
		for(int i=0; i<ret.length; ++i) {
			ret[i] = list.get(i);
		}
		return ret;
	}
}
