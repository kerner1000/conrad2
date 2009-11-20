package calhoun.analysis.crf.solver;

import java.util.List;

import calhoun.analysis.crf.FeatureList;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.util.Assert;

/** a dummy cache processor that fulfills the interface but doesn't cache.  It always retrieves values
 * by calling the {@link calhoun.analysis.crf.FeatureManager} <code>evlaute</code> functions.  <b>Does not work with 
 * NodeBoundary features.</b>
 */
public class NoCachingCacheProcessor extends CacheProcessorBasic {
	//private static final Log log = LogFactory.getLog(NoCachingCacheProcessor.class);

	boolean allPaths;
	
	/// Cached values of the sums of each feature value through the whole training set.
	double[] featureSumsLocal;
	boolean[] invalidTransitions;
	DirectFeatureList result;

	class DirectFeatureList implements FeatureList {
		FeatureEvaluation evals1;
		public int position;
		boolean valid;
		
		public DirectFeatureList() {
		}
		
		public void addFeature(int index, double val) {
			evals1.index[position] = (short) index;
			evals1.value[position++] = (float) val;
		}

		/** Returns the invalid flag. */
		public boolean isValid() {
			return valid;
		}

		/** Invalidates results. */
		public void invalidate() {
			valid = false;
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

	@Override
	public void setTrainingData(ModelManager fm, List<? extends TrainingSequence<?>> data) {
		super.setTrainingData(fm, data);
		basicInit(allPaths);
		result = new DirectFeatureList();
		
		invalidTransitions = new boolean[modelInfo.nPotentials*modelInfo.totalPositions];
		calcFeatureSums();
	}

	public boolean[] getInvalidTransitions() {
		return invalidTransitions;
	}

	@Override
	public double[] getFeatureSums() {
		return featureSumsLocal;
	}

	void calcFeatureSums() {
		featureSumsLocal = new double[modelInfo.nFeatures];
		for(int seqNum = 0; seqNum < modelInfo.nSeqs; ++seqNum) {
			TrainingSequence train = (TrainingSequence) data.get(seqNum);
			int len = train.length();
			int previousState = -1;
			for(int pos = 0; pos < len; ++pos) {
				evaluatePosition(seqNum, pos);
				int state = train.getY(pos);
				int i = 0;
				FeatureEvaluation potEval = evals[state];
				int index = potEval.index[i];
				while(index != -1) {
					float val = potEval.value[i];
					Assert.a(!Float.isNaN(val));
					featureSumsLocal[index] += val;
					index = potEval.index[++i];
				}
				
				for(int trans = 0; trans < modelInfo.nTransitions; ++trans) {
					if(modelInfo.transitionTo[trans] == state && modelInfo.transitionFrom[trans] == previousState) {
						i = 0;
						int pot = modelInfo.nStates+trans;
						potEval = evals[pot];
						index = potEval.index[i];
						while(index != -1) {
							float val = potEval.value[i];
							Assert.a(!Float.isNaN(val));
							featureSumsLocal[index] += val;
							index = potEval.index[++i];
						}
					}
				}
				previousState = state;
			}
		}
	}

	public void evaluatePosition(int seqNum, int pos) {
		InputSequence seq = data.get(seqNum);
		for(int pot=0; pot<modelInfo.nPotentials; ++pot) {
			result.evals1 = evals[pot];
			result.position = 0;
			result.valid = true;
			if(pot < modelInfo.nStates) {
				fm.evaluateNode(seq, pos, pot, result);
			}
			else {
				if(pos == 0) {
					result.evals1.index[0] = -1;
					continue;
				}
				int trans = pot - modelInfo.nStates;
				fm.evaluateEdge(seq, pos, modelInfo.transitionFrom[trans], modelInfo.transitionTo[trans], result);
			}
			if(result.isValid()) {
				result.evals1.index[result.position] = -1;
			}
			else {
				result.evals1.value[0] = Float.NaN;
				result.evals1.index[0] = Short.MIN_VALUE;
				result.evals1.index[1] = -1;
			}
		}
	}

	public void evaluateSegmentsEndingAt(int seq, int pos) {
		throw new UnsupportedOperationException();
	}

}
 