package calhoun.analysis.crf.solver.check;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.analysis.crf.solver.CacheProcessor;
import calhoun.util.Assert;

/** Extends the regular feature cache to handle explicit-length features.
 * 
 * The regular FeatureCache stays as it is.  Length based features are appended to that.
 * 
 * The length cache is like a giant array.
 * feature[seq][position][length][potential][id][val]
 * 
 * seqIndex indexes into the starts array
 * 
 * The final values are denormalized into 4 arrays:
 * length
 * potential
 * feature id
 * val
 * 
 * For every valid length and potential at least one entry in the cache will exist.  If no features are valid for a 
 * given length and potential a feature id of -1 is used.
 * 
 * The rationale behind this design is to make the cache maximally efficient for a standard case of expliciti length nodes.
 * This case is where every valid length has one or a small number of length dependent features firing and only a small set of 
 * potentials will be valid.
 * 
 * For gene calling, this has the suboptimality that intron states are duplicates of each other and will require duplication in the cache.
 * 
 * It takes 9 bytes per lookback cache entry.
 */
public class FeatureCacheLength extends FeatureCache {
	private static final Log log = LogFactory.getLog(FeatureCacheLength.class);

	/// List of states that have lookback enabled and the potentials for each of them
	public CacheProcessor.StatePotentials[] statesWithLookback;
	
	/// List of states that don't have lookback and the potentials for each of them.
	public CacheProcessor.StatePotentials[] statesWithoutLookback;

	/// Max lookback for each state
	public short[] maxStateLengths;
	public short[] minStateLengths;

	/// The starts array for the cache of length dependent features
	public int[] lengthStarts;

	/// Total size of the lookback cache
	public short maxLookback=1;
	public int nLookbackFeatures;

	public short[] lookbacks;
	public byte[] lengthPotentials;
	public short[] lengthIndexes;
	public float[] lengthVals;

	// Number of dummy features in the cache.
	int totalDummy;
	
	/// Whether a feature has been found in the current lookback
	boolean featureFound = false;
	int currentFeature;
	
	public FeatureCacheLength(ModelManager fm, List<? extends TrainingSequence<?>> data, boolean allPaths, short[] maxStateLengths, short[] minStateLengths, boolean ignoreSemiMarkovSelfTransitions) {
		this.maxStateLengths = maxStateLengths;
		this.minStateLengths = minStateLengths;
		this.ignoreSemiMarkovSelf = ignoreSemiMarkovSelfTransitions;
		init(fm, data, allPaths);
		initStatesWithLookback(maxStateLengths);

		lengthStarts = new int[totalPositions+1];

		// As in the other cache initialization, we use a first pass to estimate size and a second to fill in the values.
		computeLengthCache(true);
		nLookbackFeatures = currentFeature;
		
		lookbacks = new short[nLookbackFeatures];
		lengthPotentials = new byte[nLookbackFeatures];
		lengthIndexes = new short[nLookbackFeatures];
		lengthVals = new float[nLookbackFeatures];
		
		computeLengthCache(false);
		log.info(String.format("%d length features. %d are constraints.  Max lookback %d", nLookbackFeatures, totalDummy, maxLookback));
	}

	@Override
	protected boolean allowSelf(int state) {
		return maxStateLengths[state] > 1;
	}
	
	/** Calculates the cache for length based features.
	 * @param estimate If true, the size of the cache is calculated.  If false, the cache is actually populated. */
	void computeLengthCache(boolean estimate) {
		currentFeature = 0;
		ArrayFeatureList nodeFeatureList = new ArrayFeatureList(fm);
		ArrayFeatureList edgeFeatureList = new ArrayFeatureList(fm);
		// Iterate through all of the training sequences.  */
		for (int i = 0; i < nSeqs; ++i) {
			TrainingSequence seq = data.get(i);
			int seqLen = seq.length();

			short lastState = -1;
			short length = 0; 			
			short trainingState = -1;
			boolean segmentEnd = true;

			// Iterate through each position in the current training sequence
			for (int pos = 0; pos < seqLen; ++pos) {
				if(segmentEnd) {
					if(trainingState != -1) {
						Assert.a(length == 0 || maxStateLengths[trainingState]==1 || (length <= maxStateLengths[trainingState] && length >= minStateLengths[trainingState]), "Seq #"+i+" Pos "+pos+" Training segment "+length+" is longer than allowed length "+maxStateLengths[trainingState]);
					}
					lastState = trainingState;
					length = 0;
				}
				else {
					++length;
				}

				trainingState = (short) seq.getY(pos);
				segmentEnd = pos == seqLen-1 || trainingState != seq.getY(pos+1);
				int overallPos = seqOffsets[i] + pos;
				lengthStarts[overallPos] = currentFeature;
				// Iterate through those states that have explicit lengths
				int positionIndex = overallPos*nPotentials;
				for(CacheProcessor.StatePotentials statePotentials : statesWithLookback) {
					byte currentState = statePotentials.state;

					// First check that the state can end at this position
					boolean validExit = checkExit(positionIndex, pos, seqLen, currentState);
					if(!validExit)
						continue;

					int selfTrans = transitionIndex.getQuick(currentState, currentState) + nStates;
					Assert.a(selfTrans != -1);

					// If the state can end at this position, begin looking for possible lookback distances
					short maxStateLookback = maxStateLengths[currentState];
					// Determine the longest lookback used by any state.  Useful for iteration later on.
					maxLookback = (short) Math.max(maxLookback, maxStateLookback);
					for(short lookback = 0; lookback < maxStateLookback; ++lookback) {
						int startPos = pos - lookback;
						int lookbackPosIndex = startPos*nPotentials;
						
						// If this lookback position doesn't allow the current state, ignore this length and all longer lengths
						if (lookbackPosIndex + currentState < 0 || lookbackPosIndex + currentState >= invalidTransitions.length) {
							continue;
						}
						// Current lookback can be disallowed by node invalidation or by invalidating the self-transition.
						if(invalidTransitions[lookbackPosIndex + currentState] || (lookback>0 && invalidTransitions[lookbackPosIndex + nPotentials + selfTrans]))
							break;
						
						if(lookback+1 < minStateLengths[currentState])
							continue;
						
						// At this point start tracking if we have found a feature
						// Look for the length-dependent node features and save them.
						nodeFeatureList.clear();
						nodeFeatureList.evaluateNodeLength(seq, pos, lookback+1, currentState);
						if(!nodeFeatureList.isValid()) {
							// Exit if invalid
							continue;					
						}

						boolean computeNodeSums = segmentEnd && currentState == trainingState && lookback == length;

						// If this is a lookback that goes back to the beginning, we don't worry about edge features.
						if(startPos == 0) {
							// If we are at the beginning, we can stop here.
							loadCache(estimate, computeNodeSums, lookback, currentState, nodeFeatureList);
							if(nodeFeatureList.size() == 0) {
								dummyFeature(estimate, lookback, currentState);
							}
							
							// Don't bother looking farther back
							break;
						}

						// Otherwise, start looking for a way in to this state
						// We only want to cache the node features if there is a valid edge into this state.  nodeWritten tracks that. 
						boolean nodesWritten = false;
						
						for(byte potential : statePotentials.potentials) {
							if(invalidTransitions[lookbackPosIndex + potential]) {
								continue;
							}

							int fromState = transitionFrom[potential-nStates];
							if(fromState == currentState) {
								// Disallow self-transitions for length features
								continue;
							}
							
							// We have a valid way in and out, now track the edges.
							edgeFeatureList.clear();
							edgeFeatureList.evaluateEdgeLength(seq, pos, lookback+1, fromState, currentState);
							if(!edgeFeatureList.isValid()) {
								// Exit if invalid
								continue;
							}

							// First time you do a potential, write out the nodes (if any)
							if(nodesWritten == false) {
								loadCache(estimate, computeNodeSums, lookback, currentState, nodeFeatureList);
								nodesWritten = true;
							}
							
							// Now write the edges
							loadCache(estimate, computeNodeSums && lastState == fromState, lookback, potential, edgeFeatureList);
							
							if(edgeFeatureList.size() == 0) {
								dummyFeature(estimate, lookback, potential);
							}
						}
					}
				}
			}
		}
		lengthStarts[totalPositions] = currentFeature;
	}

	/** Checks if there is a valid transition out of a node */
	boolean checkExit(int positionIndex, int pos, int seqLen, int node) {
		if(invalidTransitions[positionIndex+node])
			return false;

		// This requires that we be in the last position or that there is a valid transition out.
		if(pos == seqLen-1) {
			return true;
		}
		
		boolean wayOut = false;
		int nextPosIndex = positionIndex + nPotentials;
		for(int transIndex=0; transIndex<nStates; ++transIndex) {
			if(transIndex == node)
				continue;
			int trans = transitionIndex.getQuick(node, transIndex);
			if(trans != -1 && !invalidTransitions[nextPosIndex + trans + nStates]) {
				wayOut = true;
				break;
			}
		}
		
		return wayOut; 
	}
	
	/** In the length cache all lookbacks are assumed to be invalid.  Only ones with cache entries are assumed valid.  When a 
	 * given lookback and edge are valid, but no features exist, a dummy feature is created.  At the start of a sequence, there are no
	 * edges, and so we just have nodes.  Therefore we may have a dummy feature for nodes corresponding to the first segment if they
	 * have no features. */
	void dummyFeature(boolean estimate, short lookback, byte potential) {
		if(estimate) {
			currentFeature += 1;
			totalDummy += 1;
		}
		else {
			//log.info(String.format("Caching dummy  feature: Lb: %d Pot: %d", lookback, potential));
			lookbacks[currentFeature] = lookback;
			lengthPotentials[currentFeature] = potential;
			lengthIndexes[currentFeature] = -1;
			currentFeature++;
		}
	}

	void loadCache(boolean estimate, boolean computeSums, short lookback, byte potential, ArrayFeatureList featureList) {
		int size = featureList.size();
		if(estimate) {
			currentFeature += size;
		}
		else {
			for(int i=0; i<size; ++i) {
				lookbacks[currentFeature] = lookback;
				lengthPotentials[currentFeature] = potential;
				short index = (short) featureList.getIndex(i);
				lengthIndexes[currentFeature] = index;
				float val = (float) featureList.getValue(i);
				lengthVals[currentFeature] = val;
				if(computeSums) {
					featureSums[index] += val;
				}
				//log.info(String.format("Caching length feature: Lb: %d Ix: %d Pot: %d Val: %f", lookback, featureList.getIndex(i), potential, lengthVals[currentFeature]));
				++currentFeature;
			}
		}
	}
	
	/** Creates an array of StatePotential objects given the maximum lookback for each state. */
	void initStatesWithLookback(short[] maxStateLengths) {
		List<CacheProcessor.StatePotentials> with = new ArrayList<CacheProcessor.StatePotentials>();
		List<CacheProcessor.StatePotentials> without = new ArrayList<CacheProcessor.StatePotentials>();
		for(byte i=0; i< maxStateLengths.length; ++i) {
			CacheProcessor.StatePotentials p = new CacheProcessor.StatePotentials();
			p.state = i;
			List<Byte> pots = new ArrayList<Byte>();
			boolean length = maxStateLengths[i] > 1;
			for(int prevState = 0; prevState < maxStateLengths.length; ++prevState) {
				int trans = transitionIndex.getQuick(prevState, i);
				if(trans != -1) {
					pots.add((byte) (trans + nStates));
				}
			}
			p.potentials = toByteArray(pots);
			if(length) {
				with.add(p);
			}
			else {
				without.add(p);
			}
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
