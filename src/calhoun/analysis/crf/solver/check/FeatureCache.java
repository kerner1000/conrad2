package calhoun.analysis.crf.solver.check;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.FeatureList;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.util.Assert;
import calhoun.util.CheckException;
import calhoun.util.ConfigException;

public class FeatureCache extends TransitionInfo {
	private static final Log log = LogFactory.getLog(FeatureCache.class);

	/// Cache feature information
	public short[] id;
	public byte[] potentialIx;
	public float[] val;
	public int longestSeq = 0; 
	public int totalPositions = 0; 
	public int totalFeatures = 0; 

	/// Values calucated as part of constant feature evaluation
	public ArrayList<Short>[] constantId;
	public ArrayList<Float>[] constantVal;
	public int cachedFeatures = 0; 
	public int numConstantFeatures = 0;
	public boolean[] invalidTransitions;
	
	/// Cached values of the sums of each feature value through the whole training set.
	public double[] featureSums;
	/// Index into the feature arrays of the first feature for each postion of each sequence.
	public int[] starts;
	/// Index into the starts array of the first position of each sequence.
	public int[] seqOffsets;

	// Basic data
	protected List<? extends TrainingSequence<?>> data;
	protected int nSeqs;

	public FeatureCache(ModelManager fm, List<? extends TrainingSequence<?>> data) {
		this(fm, data, false);
	}
	
	public FeatureCache(ModelManager fm, List<? extends TrainingSequence<?>> data, boolean allPaths) {
		init(fm, data, allPaths);
	}
	
	protected FeatureCache() {
	}
	
	protected void init(ModelManager fm, List<? extends TrainingSequence<?>> data, boolean allPaths) {
		// Constant featurallPathses and transitions
		initTrans(fm, allPaths);
		this.data = data;
		nSeqs = data.size();
		
		// Creating this object initalizes the constant features
		/* Estimate the cache size and compute the constant features */
		new ConstantFeatureList();
		totalFeatures = cachedFeatures + numConstantFeatures;
		new CachedFeatureList();
	}

	class CachedFeatureList implements FeatureList {
		// State information
		/// Current index into the cache
		private int currentFeature;
		/// The current constant feature we are checking for when evaluating the cache
		private int constantFeature;
		/// Index of the potential we are evaluating (potential is a node or edge)
		private int potentialIndex;
		/// The total number of positions in all sequences before the current one.
		private int previousPositions = 0;
		/// The previous state for the current potential being evaluated. -1 for node potentials
		private byte prevState;
		/// The current state for the current potential being evaluated.
		private byte state;
		// / True if the current potential is in the training data and so we should sum it.
		private boolean computeSum;
	
		public CachedFeatureList() {
			// Allocate cache
			featureSums = new double[fm.getNumFeatures()];

			seqOffsets = new int[nSeqs+1];
			starts = new int[totalPositions + 1];
			starts[totalPositions] = totalFeatures;
	
			id = new short[totalFeatures];
			potentialIx = new byte[totalFeatures];
			val = new float[totalFeatures];
	
			// Cache the constant features
			currentFeature = 0;
			for (short potential : orderedPotentials) {
				for (int j = 0; j < constantId[potential].size(); ++j) {
					potentialIx[currentFeature] = (byte) potential;
					id[currentFeature] = constantId[potential].get(j);
					val[currentFeature] = constantVal[potential].get(j);
					++currentFeature;
				}
			}
	
			// Loop through to cache the features
			seqOffsets[0] = 0;
			for (int i = 0; i < nSeqs; ++i) {
				TrainingSequence seq = data.get(i);
				int seqLen = seq.length();
				longestSeq = Math.max(seqLen, longestSeq);

				seqOffsets[i + 1] = seqOffsets[i] + seqLen;
				for (int pos = 0; pos < seqLen; ++pos) {
					cachePosition(seq, i, pos);
				}
				previousPositions += seq.length();
			}
	
			Assert.a(currentFeature == totalFeatures, "Current features: ", currentFeature, " Total Features: ", totalFeatures);
			log.info(String.format("Cached %d positions in %d sequences. %d features.  %d constant features.", seqOffsets[nSeqs], nSeqs, cachedFeatures,
					numConstantFeatures));
		}
		
		void cachePosition(TrainingSequence seq, int i, int pos) {
			starts[seqOffsets[i]+pos]=currentFeature;
	
			byte trainingState = (byte) seq.getY(pos);
			byte trainingPrev = (byte) (pos == 0 ? -1 : seq.getY(pos-1));
	
			int invalidIndex = (previousPositions + pos)*nPotentials;
			for(short pot : orderedPotentials) {
				constantFeature = 0;
				potentialIndex = pot;
				if(pot < nStates) {
					// Evaluate nodes
					computeSum = pot == trainingState;

					state = (byte) pot;
					prevState = (byte) -1;
					
					if(invalidTransitions[invalidIndex+state]) {
						if(computeSum) {
							throw new ConfigException("Seq: "+i+" Pos: "+pos+". Illegal state in training data: "+fm.getStateName(trainingState));
						}
						continue;
					}
	

					fm.evaluateNode(seq, pos, state, this);
				}
				else {
					if(pos == 0)
						continue;
					int trainingTransition = transitionIndex.getQuick(trainingPrev, trainingState); 
					// Verify that this transition is a valid one in the model (can still be invalid through a constraint).
					if(trainingTransition == -1) {
						throw new ConfigException("Seq: "+i+" Pos: "+pos+". Transition in training data that is disallowed in the model: "+fm.getStateName(trainingPrev)+"-"+fm.getStateName(trainingState));
					}

					int trans = pot - nStates;
					computeSum = trans == trainingTransition;
					prevState = (byte) transitionFrom[trans];
					state = (byte) transitionTo[trans];
					
					if(invalidTransitions[invalidIndex+pot]) {
						if(computeSum) {
							throw new ConfigException("Seq: "+i+" Pos: "+pos+". Transition in training data violated a constraint: "+fm.getStateName(trainingPrev)+"-"+fm.getStateName(trainingState));
						}
						continue;
					}
			

					fm.evaluateEdge(seq, pos, prevState, state, this);
				}
				
				int expectedSize = constantId[potentialIndex].size();
				Assert.a(constantFeature == expectedSize, "Pos: ",pos," had ",constantFeature," instead of ",expectedSize," constant features.");
			}
		}
		
		public void addFeature(int index, double doubleVal) {
			if (computeSum) {
				featureSums[index] += doubleVal;
			}
			float value = (float) doubleVal;
			// Check if this is the next constant feature we expect
			if (constantId[potentialIndex].size() > constantFeature) {
				// We expect a constant feature. See if this is it.
				if (constantId[potentialIndex].get(constantFeature) == index && constantVal[potentialIndex].get(constantFeature) == value) {
					++constantFeature;
					return;
				}
			}
			potentialIx[currentFeature] = (byte) potentialIndex;
			id[currentFeature] = (short) index;
			val[currentFeature] = value;
			++currentFeature;
		}
	
		public void invalidate() {
			throw new CheckException("A valid transition became an invalid transition during cache population.  Potential: "+potentialIndex);
		}
	
		public boolean isValid() {
			return true;
		}
	}

	/** This class exists to estimate the size of the feature cache required and determine the constant features. 
	 * It also creates a cache of the disallowed transitions. */ 
	class ConstantFeatureList implements FeatureList {
		/// True if the current potential is invalidated.
		private boolean invalid;
		/// Index of the potential we are evaluating (potential is a node or edge)
		private int potentialIndex;
		/// Set to true once a potential has been initialized and a set of constant features established 
		private boolean[] initialized;
		/// List of the constant features which were not found at the current position
		private ArrayList<Integer> constantFeaturesToRemove = new ArrayList<Integer>();
		/// Index into the constant feature list for the current potential of the next expected constant feature
		private int expectedConstantIndex;
		/// The number of features which need to be cached for the current potential 
		private int featuresThisEval;
		/// The total number of positions in all sequences before the current one.
		private int previousPositions = 0;
		/// The array that tracks for each potential the total number of times it is invalid.
		private int invalidCount[];

		/** Estimate the cache size and compute the constant features */
		public ConstantFeatureList() {
			int nPotentials = nStates + nTransitions;
			constantId = new ArrayList[nPotentials];
			constantVal = new ArrayList[nPotentials];
			invalidCount = new int[nPotentials];
			initialized = new boolean[nPotentials];

			for(int i=0; i<nPotentials; ++i) {
				constantId[i] = new ArrayList<Short>();
				constantVal[i] = new ArrayList<Float>();
			}

			int totalLength = 0;
			for(int i=0; i<nSeqs; ++i) {
				TrainingSequence seq = data.get(i);
				totalLength += seq.length();
			}
			
			invalidTransitions = new boolean[totalLength*nPotentials];
			
			// Iterate through sequences
			for(int i=0; i<nSeqs; ++i) {
				TrainingSequence seq = data.get(i);
				int seqLen = seq.length();

				// Iterate through positions
				for(int j=0; j<seqLen; ++j) {
					for(short pot : orderedPotentials) {
						// For the first position, no edge features
						if(j != 0 || pot < nStates) {
							checkPotential(seq, i, j, pot);
						}
					}
				}
				previousPositions += seq.length();
			}
			totalPositions = previousPositions;
		}

		public void addFeature(int index, double val) {
			if(invalid) {
				return;
			}
			float value = (float) val;
			
			if(!initialized[potentialIndex]) {
				// Now add the feature.
				constantId[potentialIndex].add((short) index);
				constantVal[potentialIndex].add((float) value);
			}
			else {
				// There are three possible scenarios
				// 1. This is the correct constant feature. 
				// 2. This is a feature that needs to be cached 
				// 2. We missed a constant feature
				int numConstants = constantId[potentialIndex].size();
				if(expectedConstantIndex >= numConstants) {
					// Must be cached, no constant features left
					//log.info("No constant features left");
					featuresThisEval++;
				}
				else {
					// We still have constant features to look for
					//log.info("Checking constant feature on potential "+potentialIndex+" "+constantId[potentialIndex].get(expectedConstantIndex)+" "+constantVal[potentialIndex].get(expectedConstantIndex));
					if(constantId[potentialIndex].get(expectedConstantIndex) == index && constantVal[potentialIndex].get(expectedConstantIndex) == value) {
						// This is the constant feature we expected
						//log.info("Found expected constant feature "+index+" "+value);
						expectedConstantIndex++;
					}
					else {
						// Look through the rest to determine if some constant features are not present or if this should be cached.
						boolean found = false;
						int i = expectedConstantIndex+1;
						for(; i < numConstants; ++i) {
							//log.info("Checking constant feature on potential "+potentialIndex+" "+constantId[potentialIndex].get(i)+" "+constantVal[potentialIndex].get(i));
							if(constantId[potentialIndex].get(i) == index && constantVal[potentialIndex].get(i) == value) {
								// This is a later constant feature
								found = true;
								break;
							}
						}
						if(found) {
							// Tag constant features for removal.
							for(int toRemove = expectedConstantIndex; toRemove < i; ++toRemove) {
								//log.info("Add to remove list "+toRemove);
								constantFeaturesToRemove.add(toRemove);
							}
							expectedConstantIndex = i+1;
						}
						else {
							// Feature was not found.  It should be cached. 
							featuresThisEval++;
							//log.info("New feature on potential "+potentialIndex+" "+index+" "+value);
						}
					}
				}
			}
		}
	
		public void invalidate() {
			invalid = true;
		}
	
		public boolean isValid() {
			return !invalid;
		}

		void checkPotential(TrainingSequence seq, int seqNum, int pos, int potential) {
			boolean node = potential < nStates;

			// Initialize the position
			invalid = false;
			featuresThisEval = 0;
			expectedConstantIndex = 0;
			constantFeaturesToRemove.clear();
			if(node) {
				potentialIndex = potential;
				fm.evaluateNode(seq, pos, potential, this);
			}
			else {
				potentialIndex = potential;
				fm.evaluateEdge(seq, pos, transitionFrom[potential - nStates], transitionTo[potential - nStates], this);
			}
			
			// Handle case where the potential was invalid
			if(invalid) {
				// First check that this doesn't occur in our training data
				if(node) {
					Assert.a(seq.getY(pos) != potential, "Seq: ",seqNum, " Pos: ", pos, ". Invalid state in training data: ", potential);
				}
				else {
					Assert.a(seq.getY(pos-1) != transitionFrom[potential - nStates] || seq.getY(pos) != transitionTo[potential - nStates], "Seq: ",seqNum, " Pos: ", pos, ". Invalid transition in training data: ", transitionFrom[potential - nStates], "-", transitionTo[potential - nStates]);
				}
				
				++invalidCount[potentialIndex];
				
				// No features to evaluate
				featuresThisEval = 0;

				// Cache the transition matrix
				invalidTransitions[(previousPositions + pos)*nPotentials + potentialIndex] = true;
				
				if(!initialized[potentialIndex]) {
					// If we started getting constant features, throw them away.
					constantId[potentialIndex].clear();
					constantVal[potentialIndex].clear();
				}
			}
			else {
				if(!initialized[potentialIndex]) {
					numConstantFeatures += constantId[potentialIndex].size();
					log.debug(String.format("Seq: %d Pos:%d, %d initial constant features for potential %d", seqNum, pos, constantId[potentialIndex].size(), potentialIndex));
					initialized[potentialIndex] = true;
				}
				else {
					// Add any constant features that are still present at the end.
					int constantSize = constantId[potentialIndex].size();
					if(expectedConstantIndex != constantSize) {
						Assert.a(expectedConstantIndex < constantSize);
						for(int i = expectedConstantIndex; i < constantSize; ++i) {
							//log.info("Add to end of remove list "+i);
							constantFeaturesToRemove.add(i);
						}
					}
					
					if(constantFeaturesToRemove.size() > 0) {
						numConstantFeatures -= constantFeaturesToRemove.size();
						// However many constant features we are short, we need to add that to the total feature count
						// Num of features to remove times the number of bases we have examined so far where the potential is valid
						// Num of bases examined with a valid potential = 
						// pos-1 bases from this sequence (-2 for edge potentials)
						// previousPositions bases from previous sequences (-i for edge potentials)
						// minus the number of invalid positions
						int newCachedFeatures = cachedFeatures + (constantFeaturesToRemove.size()) *(previousPositions+pos-(node ? 0 : seqNum+1)-invalidCount[potentialIndex]); 
						log.debug(String.format("Seq: %d Pos: %d Potential: %d Removed %d constant features.  Cached features goes from %d to %d.", seqNum, pos, potentialIndex, constantFeaturesToRemove.size(), cachedFeatures, newCachedFeatures));
						cachedFeatures = newCachedFeatures;
			
						int count = 0;
						for(Integer toRemove : constantFeaturesToRemove) {
							// This is tricky.  remove(int) must be called, not remove(Integer)
							int index = ((int) toRemove) - count;
							//log.info("Removing "+toRemove+" "+index+" "+count+" "+constantId[potentialIndex].size());
							constantId[potentialIndex].remove(index);
							constantVal[potentialIndex].remove(index);
							++count;
						}
						Assert.a(constantSize - constantFeaturesToRemove.size() == constantId[potentialIndex].size(), String.format("Removed %d, went from %d to %d", constantFeaturesToRemove.size(), constantSize, constantId[potentialIndex].size()));
					}
				}
			}		
			cachedFeatures += featuresThisEval;
		}
	}
}
