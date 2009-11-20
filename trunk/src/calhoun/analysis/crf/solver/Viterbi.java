package calhoun.analysis.crf.solver;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.CRFInference;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.analysis.crf.solver.check.FeatureCalculator;
import calhoun.analysis.crf.solver.check.TransitionInfo;
import calhoun.util.ColtUtil;
import calhoun.util.DenseBooleanMatrix2D;

public class Viterbi implements CRFInference {
	Log log = LogFactory.getLog(Viterbi.class);
	boolean debug = false; //log.isInfoEnabled();
	
	private double[] bestScore;
	private int[] backPointers;
	private boolean allPaths;

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

	public InferenceResult predict(ModelManager fm, InputSequence<?> seq, double[] lambda) {
		int numStates = fm.getNumStates();
		int len = seq.length();
		
		DenseBooleanMatrix2D transitions = fm.getLegalTransitions();
		if(allPaths || transitions == null) {
			transitions = new DenseBooleanMatrix2D(numStates, numStates);
			transitions.assign(true);
		}
		TransitionInfo t = new TransitionInfo(fm, false);
		FeatureCalculator calc = new FeatureCalculator(fm, lambda, t);

		bestScore = new double[numStates * len];
		backPointers = new int[numStates * len];
		for (int pos = 0; pos < len; pos++) {
			int posIndex = pos*numStates;
			for(int state = 0; state<numStates; ++state) {
				int index = posIndex + state;
				double nodeVal = calc.calcNodeValue(seq, pos, state);
				if(pos == 0) {
					//log.debug(String.format("Pos: %d State: %d Total: %.2f", pos, state, nodeVal));
					bestScore[index] = nodeVal;
				}
				else {
					double max = Double.NEGATIVE_INFINITY;
					int prevState = -1;
					for(int k=0; k<numStates; ++k) {
						if(!transitions.getQuick(k, state)) {
							continue;
						}
						double previous = bestScore[posIndex-numStates+k];
						double edge = calc.calcEdgeValue(seq, pos, k, state);
						double current = previous + edge + nodeVal;
						if(debug)
							log.debug(String.format("Pos: %d Trans: %d-%d %.2f (Prev: %.2f + Edge: %.2f + Node: %.2f)", pos, k, state, current, previous, edge, nodeVal));
						if(current > max) {
							//log.info("Selected Pos: "+pos+" Edge: "+k+"-"+state+" "+current+" vs. "+max+" Prev. "+previous+" Node "+nodeVal + " Calc "+calc.calcEdgeValue(seq, pos, k, state));
							max = current;
							prevState = k;
						}
					}
					//Assert.a(prevState != -1, "No legal transitions found to state ", state, ".  Pos ", pos);
					// With constraints it is legal to have states that are disallowed at a given position.
					// Just put a -infinity in there for the score and a -1 for the backpointer.
					bestScore[index] = max;
					backPointers[index] = prevState;
					//log.info(String.format("Pos: %d Prev: %d State: %d Total: "+max, pos, prevState, state));
				}
			}
		}
		
		//log.info(ColtUtil.format(bestScore));
		//log.info(backPointers);
		
		int[] ret = new int[len];
		ret[len-1] = ColtUtil.maxInColumn(bestScore, numStates, len-1);
		for(int i = len-1; i>0; --i) {
			ret[i-1] = backPointers[numStates*i + ret[i]];
		}
		InferenceResult inferenceResult = new InferenceResult();
		inferenceResult.hiddenStates = ret;
		inferenceResult.bestScores = new double[numStates];
		System.arraycopy(bestScore, numStates*(len-1), inferenceResult.bestScores, 0, numStates );
		return inferenceResult;
	}

	public int[] getBackPointers() {
		return backPointers;
	}

	public void setBackPointers(int[] backPointers) {
		this.backPointers = backPointers;
	}

	public double[] getBestScore() {
		return bestScore;
	}

	public void setBestScore(double[] bestScore) {
		this.bestScore = bestScore;
	}

}
