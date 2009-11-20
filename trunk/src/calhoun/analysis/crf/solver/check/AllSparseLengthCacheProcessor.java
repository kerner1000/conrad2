package calhoun.analysis.crf.solver.check;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.SemiMarkovSetup;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.analysis.crf.solver.CacheProcessor;
import calhoun.analysis.crf.solver.CacheProcessorBasic;
import calhoun.util.Assert;
import calhoun.util.DenseBooleanMatrix2D;
import calhoun.util.DenseIntMatrix2D;

public class AllSparseLengthCacheProcessor extends CacheProcessorBasic {
	@SuppressWarnings("unused")
	private static final Log log = LogFactory.getLog(AllSparseLengthCacheProcessor.class);

	int lookbackArraySize = -1;
	int lookbackArrayFeatureSize = -1;
	
	public boolean allPaths;
	short[] minStateLengths;
	boolean ignoreSemiMarkovSelfTransitions;
	boolean includeExplicitLengthEdges;
	
	// Non-length stuff
	private boolean[] invalidTransitions;
	private short[] id;
	private byte[] potentialIx;
	private float[] val;
	
	private int[] starts;
	//private double[] featureSums;

	// Length stuff
	private int[] lengthStarts;
	private short[] lengthLookbacks;
	private byte[] lengthPotentials;
	private short[] lengthIndexes;
	private float[] lengthVals;
	
	// State stuff
	private boolean constAtStartPos;
	private int[] nonconstantStarts;

	@Override
	public void setTrainingData(ModelManager fm, List<? extends TrainingSequence<?>> data) {
		super.setTrainingData(fm, data);
		initSequenceInfo();
		modelInfo.maxStateLengths = maxStateLengths;
		initTransitions(allPaths);
		evals = FeatureEvaluation.create(modelInfo.nPotentials, Math.max(5, modelInfo.nFeatures));
		
		// Create local references to cache values
		FeatureCache cache;
		if(maxStateLengths == null) {
			cache = new FeatureCache(fm, data, allPaths);
			modelInfo.maxLookback = 1;
		}
		else {
			FeatureCacheLength lenCache = new FeatureCacheLength(fm, data, allPaths, maxStateLengths, minStateLengths, ignoreSemiMarkovSelfTransitions);
			modelInfo.maxLookback = lenCache.maxLookback;
			lengthLookbacks = lenCache.lookbacks;
			lengthStarts = lenCache.lengthStarts;
			lengthPotentials = lenCache.lengthPotentials;
			lengthIndexes = lenCache.lengthIndexes;
			lengthVals = lenCache.lengthVals;
			modelInfo.initStatesWithLookback(maxStateLengths);
			
			if(lookbackArraySize == -1)
				lookbackArraySize = modelInfo.maxLookback+2;
			if(lookbackArrayFeatureSize == -1)
				lookbackArrayFeatureSize = Math.max(5, modelInfo.nFeatures);
			lengthEvals = LengthFeatureEvaluation.create(modelInfo.statesWithLookback, lookbackArraySize, lookbackArrayFeatureSize);
			cache = lenCache;
		}

		invalidTransitions = cache.invalidTransitions;
		id = cache.id;
		potentialIx = cache.potentialIx;
		val = cache.val;

		featureSums = cache.featureSums;
		starts = cache.starts;
		
		evaluateConstantFeatures(false);
	}

	@Override
	public void initTransitions(boolean allPaths) {
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
			Assert.a(ignoreSemiMarkovSelfTransitions || !isSemiMarkovState(i) || !transitions.getQuick(i, i), "Self transitions are not currently allowed in the semi-Markov model.", i);
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
					modelInfo.transitionIndex.setQuick(j, i, count);
					if(i == j)
						modelInfo.selfTransitions[i] = count;
					modelInfo.transitionFrom[count] = j;
					modelInfo.transitionTo[count] = i;
					++count;
				}
			}
		}
		Assert.a(count == modelInfo.nTransitions);
	}
	
	public boolean[] getInvalidTransitions() {
		return invalidTransitions;
	}

	@Override
	public double[] getFeatureSums() {
		return featureSums;
	}
	
	public void evaluatePosition(int seq, int pos) {
		if((pos == 0) != constAtStartPos) {
			evaluateConstantFeatures(pos == 0); 
		}

		int overallPos = modelInfo.seqOffsets[seq]+pos;
		int current = starts[overallPos];
		int cacheStop = starts[overallPos+1];
		byte cachedPotential = -1;
		float cachedVal = Float.NaN;
		short cachedId = -1;
		if(current < cacheStop) {
			cachedPotential = potentialIx[current];
			cachedVal = val[current];
			cachedId = id[current];
			++current;
		}
		for(int currentPotential : modelInfo.orderedPotentials) {
			int currentStart = nonconstantStarts[currentPotential];
			FeatureEvaluation potEval = evals[currentPotential];
			while(cachedPotential == currentPotential) {
				potEval.index[currentStart] = cachedId;
				potEval.value[currentStart++] = cachedVal;
				if(current < cacheStop) {
					cachedPotential = potentialIx[current];
					cachedVal = val[current];
					cachedId = id[current];
					++current;
				}
				else {
					break;
				}
			}
			potEval.index[currentStart] = -1;
		}
	}

	/**
	 * Used to iterate through the cache of length based features
	 * different lookback sizes and potentials. Used for updates in both the alpha and beta vectors.
	 * 
	 * Because the cache is so condensed, this function is very tricky.  The whole idea is that we rely on the cache being very well ordered
	 */
	public void evaluateSegmentsEndingAt(int seq, int pos) {
		int overallPosition = modelInfo.seqOffsets[seq]+pos;
		int lengthCacheStart = lengthStarts[overallPosition];
		int lengthCacheStop = lengthStarts[overallPosition + 1];

		// Get the next cache entry
		short cachedLookback = -1;
		byte cachedPotential = -1;
		if(lengthCacheStart < lengthCacheStop) {
			cachedLookback = lengthLookbacks[lengthCacheStart];
			cachedPotential = lengthPotentials[lengthCacheStart];
		}
		
		CacheProcessor.StatePotentials[] statesWithLookback = modelInfo.statesWithLookback;
		int nSemiMarkovStates = statesWithLookback.length;
		for(int stateIx=0; stateIx < nSemiMarkovStates; ++stateIx) {
			CacheProcessor.StatePotentials statePotentials = statesWithLookback[stateIx];
			LengthFeatureEvaluation[] lookbackEvals = lengthEvals[stateIx];

			short currentLookback = -1;
			int lookbackIndex = -1;
			LengthFeatureEvaluation featureEval = null;
			int node = statePotentials.state;
			int featureEvalIndex = 0;
			
			// Look for node features for this state
			while(cachedLookback != -1) {
				// First check that we are still in the right state
				int cachedState = cachedPotential;
				if(cachedPotential > modelInfo.nStates) {
					int trans = cachedPotential - modelInfo.nStates;
					cachedState = modelInfo.transitionTo[trans];
				}
				if(cachedState != node) {
					break;
				}
				
				// We are still in the same state.  Check that we are in the same lookback
				if(currentLookback == -1) {
					// Check if this is the first lookback
					currentLookback = cachedLookback;
					lookbackIndex = 0;
					featureEval = lookbackEvals[lookbackIndex];
					featureEval.lookback = currentLookback;
					featureEvalIndex = 0;
				}
				else if(currentLookback != cachedLookback) {
					// New lookback length
					// Close out the old eval
					featureEval.nodeEval.index[featureEvalIndex] = -1;
					// Set up the new one
					currentLookback = cachedLookback;
					featureEval = lookbackEvals[++lookbackIndex];
					featureEval.lookback = currentLookback;
					featureEvalIndex = 0;
				}

				// At this point we either have a node or an edge.
				if(cachedPotential > modelInfo.nStates) {
					// If we have an edge feature, assert that it is just a constraint (we don't have explicit length edge values yet). 
					Assert.a(lengthIndexes[lengthCacheStart] == -1, "ExplicitLength edge features are not supported");

					// If we don't have a node feature, add a dummy.
					if(featureEvalIndex == 0) {
						featureEval.nodeEval.index[featureEvalIndex++] = -1;
					}
				}
				else {
					// Always add in node features
					featureEval.nodeEval.index[featureEvalIndex] = lengthIndexes[lengthCacheStart];
					featureEval.nodeEval.value[featureEvalIndex] = lengthVals[lengthCacheStart];
					++featureEvalIndex;
				}
				
				// Go to the next cache element
				if(++lengthCacheStart < lengthCacheStop) {
					cachedLookback = lengthLookbacks[lengthCacheStart];
					cachedPotential = lengthPotentials[lengthCacheStart];
				}
				else {
					cachedLookback = -1;
				}
			}
			
			// Now close out the last one
			if(currentLookback != -1) {
				featureEval.nodeEval.index[featureEvalIndex] = -1;
			}
			lookbackEvals[++lookbackIndex].lookback = -1;
		}
	}
	

	void evaluateConstantFeatures(boolean atStartPos) {
		// A the start position we don't include edge features.
		constAtStartPos = atStartPos;

		// Reinitialize
		nonconstantStarts = new int[modelInfo.nPotentials];
		for(FeatureEvaluation indexArray: evals) {
			indexArray.index[0] = -1;
		}
		
		int current = 0;
		byte currentPotential = 0;
		int currentStart = 0;
		while(current < starts[0]) {
			byte cachedPotential = potentialIx[current];
			if(!atStartPos || cachedPotential < modelInfo.nStates) {
				//log.warn("Current potential: "+currentPotential+" index: "+currentPotential+" Start:"+currentStart);
				if(cachedPotential != currentPotential) {
					evals[currentPotential].index[currentStart] = -1;
					nonconstantStarts[currentPotential]=currentStart;
					currentStart = 0;
					currentPotential = cachedPotential;
				}
				FeatureEvaluation potEval = evals[currentPotential];
				potEval.index[currentStart] = id[current];
				potEval.value[currentStart++] = val[current];
			}
			++current;
		}
		evals[currentPotential].index[currentStart] = -1;
		nonconstantStarts[currentPotential]=currentStart;
	}

	public boolean isAllPaths() {
		return allPaths;
	}

	public void setAllPaths(boolean allPaths) {
		this.allPaths = allPaths;
	}

	public void setSemiMarkovSetup(SemiMarkovSetup setup) {
		maxStateLengths = setup.getMaxLengths();
		minStateLengths = setup.getMinLengths();
		ignoreSemiMarkovSelfTransitions = setup.isIgnoreSemiMarkovSelfTransitions();
	}

	public int getLookbackArrayFeatureSize() {
		return lookbackArrayFeatureSize;
	}

	public void setLookbackArrayFeatureSize(int lookbackArrayFeatureSize) {
		this.lookbackArrayFeatureSize = lookbackArrayFeatureSize;
	}

	public int getLookbackArraySize() {
		return lookbackArraySize;
	}

	public void setLookbackArraySize(int lookbackArraySize) {
		this.lookbackArraySize = lookbackArraySize;
	}
}
 