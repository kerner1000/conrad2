package calhoun.analysis.crf.solver;

import java.io.BufferedWriter;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.util.Assert;
import calhoun.util.DenseBooleanMatrix2D;
import calhoun.util.DenseIntMatrix2D;
import calhoun.util.FileUtil;

/** basic functionality common to most cache processors.  */
public abstract class CacheProcessorBasic implements CacheProcessor {
	private static final Log log = LogFactory.getLog(CacheProcessorBasic.class);
	
	String trainingFile = null;
	BufferedWriter trainingWriter = null;
	protected short[] maxStateLengths;

	protected List<? extends TrainingSequence<?>> data;
	protected ModelManager fm;
	
	protected SolverSetup modelInfo;
	protected FeatureEvaluation[] evals;
	protected LengthFeatureEvaluation[][] lengthEvals;

	protected double[] featureSums; // sums of each feature through entire training data set; length of array is number of features in fm.
	protected double[][] seqFeatureSums; // sums of each feature through each sequence; length of array is number of features in fm.

	protected void computeFeatureSums() {
		trainingWriter = FileUtil.safeOpen(trainingFile);
		
		int numFeatures = fm.getNumFeatures();
		seqFeatureSums = new double[data.size()][numFeatures];
		featureSums = new double[numFeatures];
		double[] lastSegmentFeatureSums = new double[numFeatures];
		
		for (int seqnum=0; seqnum<data.size(); seqnum++) {
			TrainingSequence seq = data.get(seqnum);

			int prevSegmentState = -1;
			int segmentLength = 1;
			int lastStart = 0;

			for (int pos=0; pos < seq.length(); pos++) {
				int state = seq.getY(pos);

				evaluatePosition(seqnum, pos);
				sumFeatures(featureSums, seqFeatureSums[seqnum], evals[state]);
				
				if (pos>0) {
					int pot = modelInfo.nStates + modelInfo.transitionIndex.getQuick(seq.getY(pos-1),state);
					sumFeatures(featureSums, seqFeatureSums[seqnum], evals[pot]);
				}
				
				// If we are at the end of a segment, evaluate length features
				if(pos == seq.length() - 1 || state != seq.getY(pos+1)) {
					evaluateSegmentsEndingAt(seqnum, pos);
					
					// Find the corrent node
					int nodeIndex=0;
					for(; nodeIndex<modelInfo.statesWithLookback.length; ++nodeIndex) {
						if(modelInfo.statesWithLookback[nodeIndex].state == state)
							break;
					}
					if(nodeIndex != modelInfo.statesWithLookback.length) {
						// Find the correct lookback entry
						int lbIndex = 0;
						while(lengthEvals[nodeIndex][lbIndex].lookback != segmentLength-1) {
							if(lengthEvals[nodeIndex][lbIndex].lookback == -1) {
								Assert.a(false, "Lookback not listed. State: ", modelInfo.statesWithLookback[nodeIndex].state, " Seq: ", seqnum, " Pos: ", pos, " Len: ", segmentLength, " # Lookbacks: ", lbIndex);
							}
							++lbIndex;
						}
						
						sumFeatures(featureSums, seqFeatureSums[seqnum], lengthEvals[nodeIndex][lbIndex].nodeEval);
						if(prevSegmentState != -1 && lengthEvals[nodeIndex][lbIndex].edgeEvals != null) {
							throw new UnsupportedOperationException("ComputeFeatureSums doesn't handle explicit length edge evals yet.");
						}
					}
					prevSegmentState = state;
					lastStart = pos-segmentLength+1;
					segmentLength = 1;
				}
				else {
					segmentLength += 1;
				}
				
				// Write out the sums for this segment
				if(trainingWriter != null) {
					if(pos == seq.length() - 1 || state != seq.getY(pos+1)) {
						for (int i=0; i<numFeatures; i++) {
							FileUtil.safeWrite(trainingWriter, String.format("Seq: %d Seg: %d-%d State: %d Feat: %d Val: %f\n", seqnum, lastStart, pos, state, i, featureSums[i] - lastSegmentFeatureSums[i]));
							lastSegmentFeatureSums[i] = featureSums[i];
						}
					}
				}
			}
		}
		if(log.isDebugEnabled()) {
			log.debug("We just computed the feature sums on the training data.  The feature sums are (id,name,sum)");
			for (int j=0; j<numFeatures; j++) {
				log.debug("(" + j + ","+fm.getFeatureName(j) + "," + featureSums[j] + ")");
			}
		}
		FileUtil.safeClose(trainingWriter);
	}

	void sumFeatures(double[] featureSums, double[] seqFeatureSum, FeatureEvaluation eval) {
		int i=0;
		while(eval.index[i] != -1) {
			featureSums[eval.index[i]] += eval.value[i];
			seqFeatureSum[eval.index[i]] += eval.value[i];
			++i;
		}
	}
	
	public void setInputData(ModelManager fm, InputSequence<?> seq) {
		// Create a dummy set of hidden states all -1
		int[] dummyHiddenStates = new int[seq.length()];
		Arrays.fill(dummyHiddenStates, Integer.MIN_VALUE);
		setTrainingData(fm, Collections.singletonList(new TrainingSequence<Object>(seq, dummyHiddenStates)));
	}

	public void setTrainingData(ModelManager fm, List<? extends TrainingSequence<?>> data) {
		this.fm = fm;
		this.data = data;
	}

	public double[] getFeatureSums() {
		return featureSums;
	}
	
	public double[][] getSequenceFeatureSums() {
		return seqFeatureSums;
	}
	
	void basicInit(boolean allPaths) {
		initSequenceInfo();
		initTransitions(allPaths);
		evals = FeatureEvaluation.create(modelInfo.nPotentials, Math.max(5, modelInfo.nFeatures));
	}

	public SolverSetup getSolverSetup() {
		return modelInfo;
	}

	public FeatureEvaluation[] getFeatureEvaluations() {
		return evals;
	}
	
	public LengthFeatureEvaluation[][] getLengthFeatureEvaluations() {
		return lengthEvals;
	}
	
	protected void initSequenceInfo() {
		// Initial basic parameters
		modelInfo = new SolverSetup();
		modelInfo.nFeatures = fm.getNumFeatures();
		modelInfo.nStates = fm.getNumStates();

		// Information about the input data
		modelInfo.nSeqs = data.size();
		modelInfo.seqOffsets = new int[modelInfo.nSeqs+1];
		modelInfo.seqOffsets[0] = 0;
		modelInfo.longestSeq = 0;
		modelInfo.totalPositions = 0;
		for (int i = 0; i < modelInfo.nSeqs; ++i) {
			TrainingSequence seq = data.get(i);
			int seqLen = seq.length();
			modelInfo.longestSeq = Math.max(seqLen, modelInfo.longestSeq);
			modelInfo.seqOffsets[i + 1] = modelInfo.seqOffsets[i] + seqLen;
			modelInfo.totalPositions += seqLen;
		}
	}
	
	protected boolean isSemiMarkovState(int state) {
		return maxStateLengths == null ? false : maxStateLengths[state] > 1;
	}

	protected void initTransitions(boolean allPaths) {
		// Initial the transition information
		modelInfo.transitionIndex = new DenseIntMatrix2D(modelInfo.nStates, modelInfo.nStates);
		modelInfo.transitionIndex.assign(-1);
		modelInfo.selfTransitions = new int[modelInfo.nStates];
		Arrays.fill(modelInfo.selfTransitions, -1);
		DenseBooleanMatrix2D transitions = fm.getLegalTransitions();
		if(transitions == null || allPaths) {
			transitions = new DenseBooleanMatrix2D(modelInfo.nStates, modelInfo.nStates);
			transitions.assign(true);
		}

		short count = 0;
		for(short i = 0; i<modelInfo.nStates; ++i) {
			for(short j = 0; j<modelInfo.nStates; ++j) {
				if(transitions.getQuick(i, j) || (i == j && isSemiMarkovState(i)))
					count++;
			}
		}
		modelInfo.nTransitions = count;
		modelInfo.nPotentials = modelInfo.nStates + modelInfo.nTransitions;
		modelInfo.orderedPotentials = new short[modelInfo.nPotentials];
		modelInfo.transitionFrom = new short[modelInfo.nTransitions];
		modelInfo.transitionTo = new short[modelInfo.nTransitions];
		count = 0;
		int orderedCount = 0;
		for(short i = 0; i<modelInfo.nStates; ++i) {
			modelInfo.orderedPotentials[orderedCount] = i;
			orderedCount++;
			for(short j = 0; j<modelInfo.nStates; ++j) {
				if(transitions.getQuick(j, i) || (i == j && isSemiMarkovState(i))) {
					modelInfo.orderedPotentials[orderedCount] = (short) (modelInfo.nStates+count);
					orderedCount++;
					if(i == j)
						modelInfo.selfTransitions[i] = count;
					modelInfo.transitionIndex.setQuick(j, i, count);
					modelInfo.transitionFrom[count] = j;
					modelInfo.transitionTo[count] = i;
					++count;
				}
			}
		}
		Assert.a(count == modelInfo.nTransitions);
	}

	public String getTrainingFile() {
		return trainingFile;
	}

	public void setTrainingFile(String trainingFile) {
		this.trainingFile = trainingFile;
	}

	/**
	 * @return the data
	 */
	public List<? extends TrainingSequence<?>> getData() {
		return data;
	}
}
