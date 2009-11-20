package calhoun.analysis.crf.solver;

import java.util.ArrayList;
import java.util.List;

import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.util.Assert;
import calhoun.util.DenseIntMatrix2D;

/** interface to implementations of feature caches.  The <code>CacheProcessor</code> is the 
 * interface between most of the objective functions and the actual feature managers.
 */  
public interface CacheProcessor {
	
	/** This class holds the feature evaluations for a given position, or position/length combination.
	 * There are two arrays, one to hold the feature indicies, and the other to hold the values.  Each
	 * array is two-dimensional.  The first dimension is fixed length and contains one entry for each potential 
	 * (node or edge) in the model.  The potentials are in model order, with each node followed by all transitions into that node.
	 * The second dimension is variable length and contains a information on feature evaluations for that potential.
	 * In the second dimension, each list of evaluations is terminated by a -1 in the index column.  
	 */ 
	public static class FeatureEvaluation {
		FeatureEvaluation(int nFeatures) {
			index = new short[nFeatures];
			value = new float[nFeatures];
			index[0] = -1;
		}
		public short[] index;
		public float[] value;

		public static FeatureEvaluation[] create(int nPotentials, int featureSize) {
			FeatureEvaluation[] potential = new FeatureEvaluation[nPotentials];
			for(int i=0; i<potential.length; ++i) {
				potential[i] = new FeatureEvaluation(featureSize);
			}
			return potential;
		}
	}
	
	public static class LengthFeatureEvaluation {
		public LengthFeatureEvaluation(int nFeatures) {
			nodeEval = new FeatureEvaluation(nFeatures);
		}
		
		public short lookback = -1;
		public FeatureEvaluation nodeEval;
		public FeatureEvaluation[] edgeEvals;
		
		public static LengthFeatureEvaluation[][] create(StatePotentials[] statePotentials, int nLookbacks, int nFeatures) {
			LengthFeatureEvaluation[][] lookbacks = new LengthFeatureEvaluation[statePotentials.length][nLookbacks];
			for(int i=0; i<statePotentials.length; ++i) {
				for(int j=0; j<nLookbacks; ++j) {
					lookbacks[i][j] = new LengthFeatureEvaluation(nFeatures);
				}
			}
			return lookbacks;
		}
	}
	
	public static class SolverSetup {
		public int nFeatures;
		public int nStates;
		public int nPotentials;
		public int nTransitions;
		
		/// Index into the starts array of the first position of each sequence.
		public int[] seqOffsets;

		/// Number of sequences in the training data set
		public int nSeqs; 

		/// The length of the longest single sequence in the data set.
		public int longestSeq; 

		public short[] transitionFrom;
		public short[] transitionTo;
		public short[] orderedPotentials;
		public int totalPositions;
		public int[] selfTransitions;

		public DenseIntMatrix2D transitionIndex;
		// For each array, contains a list of the possible destinate
		public byte[][] exitTransitions;
		
		public short[] maxStateLengths;

		public short maxLookback = 1;
		public StatePotentials[] statesWithLookback;
		public StatePotentials[] statesWithoutLookback;

		public void setup(ModelManager fm, List<? extends TrainingSequence<?>> data, boolean allPaths, short[] maxStateLengths2, boolean ignoreSemiMarkovSelfTransitions) {
			this.maxStateLengths = maxStateLengths2;
			initStatesWithLookback(maxStateLengths);
			
			Assert.a(statesWithLookback != null);
			Assert.a(statesWithoutLookback != null);			
		}
		
		protected boolean allowSelf(int state) {
			return maxStateLengths[state] > 1;
		}
		

		/** Creates an array of StatePotential objects given the maximum lookback for each state. */
		public void initStatesWithLookback(short[] maxStateLengths) {
			List<CacheProcessor.StatePotentials> with = new ArrayList<CacheProcessor.StatePotentials>();
			List<CacheProcessor.StatePotentials> without = new ArrayList<CacheProcessor.StatePotentials>();
			exitTransitions = new byte[maxStateLengths.length][];

			for(byte i=0; i< maxStateLengths.length; ++i) {
				// Set up the statePotentials structure
				CacheProcessor.StatePotentials p = new CacheProcessor.StatePotentials();
				p.state = i;
				List<Byte> pots = new ArrayList<Byte>();
				List<Byte> exits = new ArrayList<Byte>();
				boolean length = maxStateLengths[i] > 1;
				for(int prevState = 0; prevState < maxStateLengths.length; ++prevState) {
					// Get list of transitions into this state
					int trans = transitionIndex.getQuick(prevState, i);
					if(trans != -1) {
						pots.add((byte) (trans + nStates));
					}

					// Get list of transitions out  this state
					trans = transitionIndex.getQuick(i, prevState);
					if(trans != -1 && prevState != i) {
						exits.add((byte) (trans + nStates));
					}
				}
				p.potentials = toByteArray(pots);
				if(length) {
					with.add(p);
				}
				else {
					without.add(p);
				}
				exitTransitions[i] = toByteArray(exits);
			}
			statesWithLookback = with.toArray(new CacheProcessor.StatePotentials[with.size()]);
			statesWithoutLookback = without.toArray(new CacheProcessor.StatePotentials[without.size()]);
		}

		private byte[] toByteArray(List<Byte> list) {
			byte[] ret = new byte[list.size()];
			for(int i=0; i<ret.length; ++i) {
				ret[i] = list.get(i);
			}
			return ret;
		}
	}

	public static final class StatePotentials {
		public byte state;
		public byte[] potentials;
	}

	public void setInputData(ModelManager fm, InputSequence<?> data);
	
	public void setTrainingData(ModelManager fm, List<? extends TrainingSequence<?>> data);
	public List<? extends TrainingSequence<?>> getData();
	
	public double[] getFeatureSums();
	public double[][] getSequenceFeatureSums();
	public boolean[] getInvalidTransitions();
	
	public SolverSetup getSolverSetup();
	
	public FeatureEvaluation[] getFeatureEvaluations();
	
	public LengthFeatureEvaluation[][] getLengthFeatureEvaluations();
	
	public void evaluatePosition(int seq, int pos);

	public void evaluateSegmentsEndingAt(int seq, int pos);
}
