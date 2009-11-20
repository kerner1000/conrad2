package calhoun.analysis.crf.solver;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import calhoun.analysis.crf.CacheStrategySpec;
import calhoun.analysis.crf.CompositeFeatureManager;
import calhoun.analysis.crf.FeatureManager;
import calhoun.analysis.crf.FeatureManagerEdge;
import calhoun.analysis.crf.FeatureManagerNode;
import calhoun.analysis.crf.FeatureManagerNodeBoundaries;
import calhoun.analysis.crf.FeatureManagerNodeExplicitLength;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.SemiMarkovSetup;
import calhoun.analysis.crf.CacheStrategySpec.CacheStrategy;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.analysis.crf.solver.check.ArrayFeatureList;
import calhoun.util.Assert;
import calhoun.util.CheckException;

/** A policy based cache processor.  This is the main cache processor used.  It uses the cache
 * policies specified by the feature managers to efficent cache all feature values.
 * The distinct caching strategies recognized by this CacheProcessor are:
 * <ul>
 *  <ol> <b>COMPOSITE</b> - must drill down further to a non-composite FeatureManager to learn the strategy
 *  <ol> <b>CONSTANT</b> - Triplets (potential, featureIndex, value) that do not depend on position
 *  <ol> <b>DENSE</b> - Triplets (potential, featureIndex, lookupTable), where lookupTable is of size totalPositions and different
 *    doublets (potential, featureIndex), if they always have the same value, can share the same lookupTable to save space.
 *  <ol> <b>SPARSE</b> - Quadruplets (cumulativePosition, potential, featureIndex, value) of nonzero evaluation stored in order first by cumulativePosition,
 *    and then arbitrarily, similar to the way the previous "FeatureCache" worked but without need to sub-order by potentials in "model order". 
 *  <ol> <b>LENGTHFUNCTION</b> - Quadruplets (length, state, featureIndex, value). 
 * <ol> Otherwise, treated as SPARSE
 * </ul>
 */
public class CacheProcessorDeluxe extends CacheProcessorBasic {
	private static final Log log = LogFactory.getLog(CacheProcessorDeluxe.class);
	private final static Logger LOGGER = LoggerFactory.getLogger(CacheProcessorDeluxe.class);
		
	private boolean ignoreInvalidTrainingData;
	private boolean allPaths;
	private boolean validHiddenStates = false;
	
	public boolean[] invalidTransitions;  // of size totalLength * nPotentials
	
	private int[] numFixedEvalIndices;  // of length nPotentials
	
	private boolean        forceStrategy;
	private CacheStrategy  onlyStrategy;
	
	int lookbackArraySize = -1;
	int lookbackArrayFeatureSize = -1;
	
	///////////////////// Below holds raw information for Sparsely Cached Feature Evaluations
	private ArrayList<FeatureManager> sparseFMList;
	private ArrayList<FeatureManager> lengthFMList;
	private int[] sparseStarts;	 // of length (totalPositions + 1)
	private ArrayList<SparseEvaluation> sparseEvals;
	private class SparseEvaluation {
		short featureIndex;
		int   potential;
		float value;
	}
	// This is a global variable used only by evaluation of sparse features from cache
	private int[] currentEvalIndices;   // of length nPotentials
	
	
	//////////////////// Below holds raw information for the Constant Features
	private boolean writtenYet;
	private ArrayList<ConstantEvaluation> constEvals;
	private class ConstantEvaluation {
		short featureIndex;
		int   potential;
		float value;
		int   evalIndex;
	}
	
	
	//////////////////////////////////////// Below holds raw information for DenseFeatureTables
	private ArrayList<float[]> denseTables; // first dimension is which table; second dimension is of length modelInfo.totalPositions; faster to flip-flop?
	// ********** Above variable is a memory hog
	private ArrayList<DenseEvalIndices> denseEvalsList; // use this until finilization, at which time pack into array below with identical info
	private class DenseEvalIndices {
		public float[]           lookupTable;
		public int               evalIndex;
		public FeatureEvaluation evalTable;   // a pointer to one of the tables in denseTables, several DenseEvalIndices can point to same table; the number of possible tables to choose from is nPotentials
		public int               potential;   // evalTable happens to be equal to evals.get(potential)
		public short             featureIndex;
		
	}
	
	
	//////////////// Below holds info for DenseNodeBoundary feature caches
	private ArrayList<float[]> denseBoundaryTables;  // note that these tables need 1 extra unit of lengh for each sequence, so they are of length modelInfo.totalPositions + modelInfo.nSeqs
	ArrayList<DenseNodeBoundaryIndices>[] denseNodeBoundaryEvalsList; // have a list for each potential

	private short[] minStateLengths;

	private boolean ignoreSemiMarkovSelfTransitions;
	private class DenseNodeBoundaryIndices {
		public int featureIndex;
		public float[] lookupTable;
		public int rightPad;
		public int leftPad;
		public int potential;	
	}
	
	//////////////// Below holds info for Length feature caches
	private LengthOnlyEvaluation[][] lengthTables;
	
	private class LengthOnlyEvaluation {
		public short[] featureIndex;
		public float[] value;
	}
	
	////////////////////////////////////////////
	////// MEMBER VARIABLES ABOVE, METHODS BELOW
	////////////////////////////////////////////
	
	public CacheProcessorDeluxe() {
		super();
		log.debug("Calling Cache Processor Deluxe constructor");
		forceStrategy = false;
	}
	
	public CacheProcessorDeluxe(CacheStrategy strategy) {
		super();
		log.warn("Calling CacheProcessorDeluxe constructor, specifying a strategy");
		forceStrategy = true;
		onlyStrategy = strategy;
		log.info("Constructed a cache processor and requiring it to always use the following cache strategy:");
		switch (strategy) {
		case CONSTANT:
			log.warn("CONSTANT");
		break;
		case DENSE:
			log.warn("DENSE");
		break;
		case SPARSE:
			log.warn("SPARSE");
		break;
		default:
			Assert.a(false,"ERROR - case not dealt with yet.");
		}
	}
	
	@Override
	public void setTrainingData(ModelManager fm, List<? extends TrainingSequence<?>> data) {
		super.setTrainingData(fm, data);
		basicInit(allPaths);

		if (maxStateLengths == null) {
			log.debug("maxStateLengths was not set, setting all to length 1");
			maxStateLengths = new short[fm.getNumStates()];
			Arrays.fill(maxStateLengths, (short) 1);
		}
		
		if (minStateLengths == null) {
			log.info("minStateLengths was not set, setting all to length 1");
			minStateLengths = new short[fm.getNumStates()];
			Arrays.fill(maxStateLengths, (short) 1);
		}
		
		modelInfo.setup(fm, data, allPaths, maxStateLengths, ignoreSemiMarkovSelfTransitions);
		
		modelInfo.maxLookback = 1;
		for(int stateLen : maxStateLengths) {
			modelInfo.maxLookback = (short) Math.max(stateLen, modelInfo.maxLookback);
		}

		if(lookbackArraySize == -1)
			lookbackArraySize = modelInfo.maxLookback+2;
		if(lookbackArrayFeatureSize == -1)
			lookbackArrayFeatureSize = Math.max(5, modelInfo.nFeatures);
		lengthEvals = LengthFeatureEvaluation.create(modelInfo.statesWithLookback, lookbackArraySize, lookbackArrayFeatureSize);

		validHiddenStates = data.get(0).getY(0) >= 0; 
		if(validHiddenStates) {
			// Verify that there are no invalid transitions in the training data.
			List<TrainingSequence<?>> goodData = new ArrayList<TrainingSequence<?>>();
			for (int i=0; i<data.size(); ++i) {
				TrainingSequence seq = data.get(i);
				short length = 1; 			
				int lastState = -1;
				boolean validSequence = true;
				try {
					for (int pos=0; pos < seq.length(); pos++) {
						int trainingState = seq.getY(pos);
						boolean segmentEnd = pos == seq.length()-1 || trainingState != seq.getY(pos+1);
						if(segmentEnd) {
							Assert.a(maxStateLengths[trainingState]==1 || (length <= maxStateLengths[trainingState] && length >= minStateLengths[trainingState]), "Seq #", i, " Pos ", pos, " Training segment length ", length, " state: "+trainingState+" outside the allowed length ", minStateLengths[trainingState], "-", maxStateLengths[trainingState]);
							Assert.a(lastState == -1 || modelInfo.transitionIndex.getQuick(lastState,trainingState) != -1);
							lastState = trainingState;
							length = 1;
						}
						else {
							++length;
						}
					}
				}
				catch(CheckException ex) {
					if(ignoreInvalidTrainingData) {
						log.warn("Discarding sequence "+i+" "+ex.getMessage());
						validSequence = false;
					}
					else {
						throw ex;
					}
				}
				if(validSequence) {
					goodData.add(seq);
				}
			}
			Assert.a(goodData.size() != 0, "All training sequences were invalid.");
			int diff = data.size() - goodData.size();
			if(diff != 0) {
				log.warn("Using "+goodData.size()+" training sequences.  Discarded "+diff+" because of state length or transition problems.");
				// Redo the initialization with the new set of training data
				setTrainingData(fm, goodData);
				return;
			}
		}
		
		initializeCacheProcessor();

		updateCacheProcessor(fm);
		cacheLengthFeatureManagers();
		cacheSparseFeatureManagers();
		
		evaluateConstantFeatures();

		if(validHiddenStates) {
			checkConstraintsInTrainingData();
			computeFeatureSums();
		}
	}


	private void initializeCacheProcessor() {
		debug("nPotentials="+modelInfo.nPotentials + ", totalPositions="+modelInfo.totalPositions + ", multiplied="+modelInfo.nPotentials*modelInfo.totalPositions);
		invalidTransitions = new boolean[modelInfo.nPotentials*modelInfo.totalPositions];
		
		//featureSums = new double[fm.getNumFeatures()];
		
		numFixedEvalIndices = new int[modelInfo.nPotentials];
		currentEvalIndices = new int[modelInfo.nPotentials];
		
		denseEvalsList  = new ArrayList<DenseEvalIndices>();
		denseTables = new ArrayList<float[]>(); 
		
		//denseNodeBoundaryEvalsList = new ArrayList<DenseNodeBoundaryIndices>[modelInfo.nPotentials];
		denseBoundaryTables = new ArrayList<float[]>();
		denseNodeBoundaryEvalsList = new ArrayList[modelInfo.nPotentials];
		for (int pot=0; pot<modelInfo.nPotentials; pot++) {
			denseNodeBoundaryEvalsList[pot] = new 	ArrayList<DenseNodeBoundaryIndices>();
		}
		
		constEvals      = new ArrayList<ConstantEvaluation>();
		writtenYet = false;
		
		sparseFMList = new ArrayList<FeatureManager>();
		sparseStarts = new int[modelInfo.totalPositions + 1];
		sparseEvals     = new ArrayList<SparseEvaluation>();
		
		lengthFMList = new ArrayList<FeatureManager>();
		lengthTables = new LengthOnlyEvaluation[modelInfo.nStates][];
		for(int i= 0; i<modelInfo.nStates; ++i) {
			lengthTables[i] = new LengthOnlyEvaluation[modelInfo.maxStateLengths[i]];
			for(int j=0; j<lengthTables[i].length; ++j) {
				lengthTables[i][j] = new LengthOnlyEvaluation();
			}
		}
	}
	
	
	private void updateCacheProcessor(FeatureManager fm1) {
		log.debug("Calling updateCacheProcessor");
		CacheStrategySpec.CacheStrategy strategy;
		if (forceStrategy) {
			strategy = onlyStrategy;
		} else {
			strategy = fm1.getCacheStrategy().strategy;
		}
		switch (strategy) {
		case COMPOSITE:
			log.debug("  burrowing into a COMPOSITE feature manager");
		CompositeFeatureManager cfm = (CompositeFeatureManager) fm1;
		for (FeatureManager fmchild : cfm.getComponentFeatures()) {
			updateCacheProcessor(fmchild);
		}
		break;
		case DENSE:
			log.debug("  processing a feature manager as DENSE");
			cacheFeatureManagerAsDenseNode(fm1);
		break;
		case DENSE_NODE_BOUNDARY:
			log.debug("Processing a node as DENSE_NODE_BOUNDARY");
			cacheFeatureManagerAsDenseNodeBoundary(fm1);
		break;
		case CONSTANT:
			log.debug("  processing a feature manager as CONSTANT");
			cacheFeatureManagerAsConstant(fm1);
		break;
		case LENGTHFUNCTION:
			log.debug("  processing a feature manager as LENGTH");
			lengthFMList.add(fm1);
		break;
		case SPARSE:
		case UNSPECIFIED:
		default: 
			log.debug("  processing a feature manager as SPARSE");
			sparseFMList.add(fm1);
			if(fm1 instanceof FeatureManagerNodeExplicitLength) {
				lengthFMList.add(fm1);
			}
		break;
		}
	}
	
	
	
	InputSequence getSequence(InputSequence seq, FeatureManager argFm) {
		return argFm.getInputComponent() == null ? seq : seq.getComponent(argFm.getInputComponent());
	}
	
	private void cacheSparseFeatureManagers() {
		sparseStarts[0] = 0;
		
		log.debug("We are caching this many feature managers (some may be composite) as sparse: " + sparseFMList.size());
		log.debug("Number of feature in each of the FMs is");
		for (FeatureManager fmb : sparseFMList) {
			log.debug("Number of features is " + fmb.getNumFeatures());
			if (FeatureManagerEdge.class.isInstance(fmb))  {log.debug("  Is edge featuremanager"); }
			if (FeatureManagerNode.class.isInstance(fmb))  {log.debug("  Is node featuremanager"); }
		}
		
		ArrayFeatureList result = new ArrayFeatureList(fm);				
		for (int i=0; i<data.size(); i++) {
			TrainingSequence seq = data.get(i);
			for (int pos=0; pos<seq.length(); pos++) {
				int cumulativePosition = modelInfo.seqOffsets[i] + pos;
				
				for (int potential=0; potential< modelInfo.nPotentials; potential++) {
					result.clear();				
					
					int state;
					int prevState;
					
					
					if (potential<modelInfo.nStates) {
						state = potential;
						for (FeatureManager fmb : sparseFMList) {
							if (!FeatureManagerNode.class.isInstance(fmb)) { continue; }							
							
							FeatureManagerNode fmn = (FeatureManagerNode) fmb;
							fmn.evaluateNode(getSequence(seq, fmb),pos,state,result);
						}
					} else if (pos>0) {
						prevState = modelInfo.transitionFrom[potential - modelInfo.nStates];
						state = modelInfo.transitionTo[potential - modelInfo.nStates];
						for (FeatureManager fm1 : sparseFMList) {
							if (!FeatureManagerEdge.class.isInstance(fm1))  {continue; }
							FeatureManagerEdge fme = (FeatureManagerEdge) fm1;
							
							fme.evaluateEdge(getSequence(seq, fme),pos,prevState,state,result);				
						}
					}
					
					
					if (result.isValid()) {
						for (int j=0; j< result.size(); j++) {
							float value = (float) result.getValue(j); 
							if(value != 0.0) {
								SparseEvaluation se = new SparseEvaluation();
								se.featureIndex = (short) result.getIndex(j);
								se.potential = potential;
								se.value = value;
								sparseEvals.add(se);
							}
						}
					} else {
						int ind = cumulativePosition*modelInfo.nPotentials + potential;
						invalidTransitions[ind] = true;
					}
				}
				sparseStarts[cumulativePosition+1] = sparseEvals.size();
			}
		}
	}
	
	private void cacheFeatureManagerAsConstant(FeatureManager fm1) {
		CacheStrategySpec s = (CacheStrategySpec) fm1.getCacheStrategy();
		if (s.strategy != CacheStrategy.CONSTANT) {
			log.warn("Caching a Feature Manager as CONSTANT even though it requested a different strategy");
		}
		
		InputSequence seq = data.get(0);
		int pos = 0;
		
		ArrayFeatureList result = new ArrayFeatureList(fm);
		for (int pot=0; pot< modelInfo.nPotentials; pot++) {
			result.clear();
			
			int state;
			int prevState;
			
			if (pot<modelInfo.nStates) {
				if (!FeatureManagerNode.class.isInstance(fm1)) { continue; }
				state = pot;
				FeatureManagerNode fmn = (FeatureManagerNode) fm1;
				fmn.evaluateNode(getSequence(seq, fmn),pos,state,result);				
			} else {
				if (!FeatureManagerEdge.class.isInstance(fm1)) { continue; }
				FeatureManagerEdge fme = (FeatureManagerEdge) fm1;
				prevState = modelInfo.transitionFrom[pot - modelInfo.nStates];
				state = modelInfo.transitionTo[pot - modelInfo.nStates];
				fme.evaluateEdge(getSequence(seq, fme),pos,prevState,state,result);
			}
			
			Assert.a(result.isValid(),"Constant features should not have to encounter invalid potentials");

			for (int j=0; j<result.size(); j++) {
				ConstantEvaluation ce = new ConstantEvaluation();
				ce.featureIndex = (short) result.getIndex(j);
				ce.potential = pot;
				ce.value = (float) result.getValue(j);
				ce.evalIndex = numFixedEvalIndices[pot];
				numFixedEvalIndices[pot]++;
				constEvals.add(ce);
			}
		}
	}
	
	class FeaturePotential {
		int featureIndex;
		int potential;
		@Override
		public int hashCode() {
			return potential * modelInfo.nPotentials + featureIndex;
		}
		@Override
		public boolean equals(Object rhs) {
			FeaturePotential f2 = (FeaturePotential) rhs;
			return featureIndex == f2.featureIndex && potential == f2.potential;
		}
	}
	
	private void cacheFeatureManagerAsDenseNode(FeatureManager fm1) {
		CacheStrategySpec s = (CacheStrategySpec) fm1.getCacheStrategy();
		if (!(s.strategy == CacheStrategy.DENSE)) {
			log.warn("The evaluations are being cached using the DENSE strategy even though that was not the strategy requested by the feature.");
		}
		
		CacheStrategySpec.DenseCachingDetails details = (CacheStrategySpec.DenseCachingDetails) s.details;
		
		if (details == null) { 
			if(fm1.getNumFeatures() > 0) {
				log.debug("No details specified for "+fm1+".  Using evaluations to find the correct nodes");
			}
			
			// Iterate through all sequences and positions to find the feature/potential combinations used
			// This replaces code that used to do it based on the first position of the first sequence with length > 2.
			HashSet<FeaturePotential> usedPotentials = new HashSet();
			ArrayFeatureList result = new ArrayFeatureList(fm);
			FeaturePotential current = new FeaturePotential();
			for (InputSequence seq : data) {
				for(int pos = 0; pos < seq.length(); ++pos) {
					for (int pot=0; pot< modelInfo.nPotentials; pot++) {
						result.clear();
						if (pot<modelInfo.nStates) {
							if (!FeatureManagerNode.class.isInstance(fm1)) { continue; }
							FeatureManagerNode fmn = (FeatureManagerNode) fm1;
							fmn.evaluateNode(getSequence(seq, fmn),pos,pot,result);				
						} else {
							if (!FeatureManagerEdge.class.isInstance(fm1)) { continue; }
							FeatureManagerEdge fme = (FeatureManagerEdge) fm1;
							int prevState = modelInfo.transitionFrom[pot - modelInfo.nStates];
							int state = modelInfo.transitionTo[pot - modelInfo.nStates];
							fme.evaluateEdge(getSequence(seq, fme),pos,prevState,state,result);
						}
						for (int j=0; j<result.size(); j++) {
							current.featureIndex = result.getIndex(j); 
							current.potential = pot;
							if(!usedPotentials.contains(current)) {
								log.debug("Adding dense node - pot: "+current.potential+" Feat: "+current.featureIndex);
								usedPotentials.add(current);
								current = new FeaturePotential();
							}
						}
					}
				}
			}

			// Now that we have the list of potentials that will be used, go and create the details.
			details = new CacheStrategySpec.DenseCachingDetails();
			ArrayList<Integer> densePotList = new ArrayList<Integer>();
			ArrayList<Integer> denseTabList = new ArrayList<Integer>();			
			ArrayList<Short> denseFiList = new ArrayList<Short>();
			details.nEvals = usedPotentials.size();
			details.nTables = 0;

			for(FeaturePotential fp : usedPotentials) {
				densePotList.add(fp.potential);
				denseTabList.add(details.nTables);
				denseFiList.add((short) fp.featureIndex);
				details.nTables++;
			}		

			details.featureIndex = new short[details.nEvals];
			details.tableNum     = new int[details.nEvals];
			details.potential    = new int[details.nEvals];
			
			for (int j=0; j<details.nEvals; j++) {
				details.featureIndex[j] = denseFiList.get(j);
				details.potential[j] = densePotList.get(j);
				details.tableNum[j] = denseTabList.get(j);
			}
			details.check();
		} 

		int firstTableIndex = denseTables.size();
		details.check();
		for (int j=0; j<details.nTables; j++) {
			float[] temp = new float[modelInfo.totalPositions];
			for (int k=0; k<modelInfo.totalPositions; k++) {
				temp[k] = 0;
			}
			denseTables.add(temp);
		}
		
		for (int j=0; j<details.nEvals; j++) {
			DenseEvalIndices de = new DenseEvalIndices();
			
			int potential = details.potential[j];
			de.evalTable = evals[potential];
			de.potential = potential;
			de.evalIndex = numFixedEvalIndices[potential];
			numFixedEvalIndices[potential]++;
			de.lookupTable = denseTables.get(firstTableIndex + details.tableNum[j]);
			de.featureIndex = details.featureIndex[j];
			denseEvalsList.add(de);
		}
		
		
		ArrayFeatureList result = new ArrayFeatureList(fm);				
		for (int i=0; i<data.size(); i++) {
			TrainingSequence seq = data.get(i);
			int actualPreviousState = -1;
			for (int pos=0; pos<seq.length(); pos++) {
				int actualState = seq.getY(pos);
				
				int cumulativePosition = modelInfo.seqOffsets[i] + pos;
				
				for (int potential=0; potential< modelInfo.nPotentials; potential++) {
					result.clear();				
					
					int state;
					int prevState;
					
					
					if (potential<modelInfo.nStates) {
						if (!FeatureManagerNode.class.isInstance(fm1)) { continue; }
						state = potential;
						FeatureManagerNode fmn = (FeatureManagerNode) fm1;
						fmn.evaluateNode(getSequence(seq, fmn),pos,state,result);
						if (state==actualState) {
							for (int j=0; j< result.size(); j++) {
//								featureSums[result.getIndex(j)] += result.getValue(j); 
							}								
						}
					} else if (pos>0) {
						if (!FeatureManagerEdge.class.isInstance(fm1)) { continue; }
						FeatureManagerEdge fme = (FeatureManagerEdge) fm1;
						prevState = modelInfo.transitionFrom[potential - modelInfo.nStates];
						state = modelInfo.transitionTo[potential - modelInfo.nStates];
						fme.evaluateEdge(getSequence(seq, fme),pos,prevState,state,result);
						if ((state==actualState) && (prevState==actualPreviousState)) {
							for (int j=0; j< result.size(); j++) {
//								featureSums[result.getIndex(j)] += result.getValue(j); 
							}								
						}
					}
					
					
					if (result.isValid()) {
						
						for (int j=0; j<result.size(); j++) {
							for (int ev=0; ev<details.nEvals; ev++) {
								DenseEvalIndices de = denseEvalsList.get(denseEvalsList.size()-1-ev);
								if ((potential==de.potential) && (result.getIndex(j)==de.featureIndex)) {
									de.lookupTable[modelInfo.seqOffsets[i] + pos] = (float) result.getValue(j);
									break;
								}
							}
						}
						
						
					} else {
						int ind = cumulativePosition*modelInfo.nPotentials + potential;
						invalidTransitions[ind] = true;
					}
				}
				actualPreviousState = actualState;
			}
		}
		
		
	}
	
	
	private void cacheFeatureManagerAsDenseNodeBoundary(FeatureManager fm1) {
		// NOTE: This does not calculate feature sums for the semi-Markov features.  I think it would be best
		// to have a separate function do this which calls the CacheProcessor but which is separate and not interwoven with each part of CacheProcessor.
		
		CacheStrategySpec s = (CacheStrategySpec) fm1.getCacheStrategy();
		Assert.a(s.strategy == CacheStrategy.DENSE_NODE_BOUNDARY,"ERROR: Cannot cache using thie DENSE_NODE_BOUNDARY strategy unless specifically requested by the feature.");
		
		CacheStrategySpec.DenseBoundaryCachingDetails details = (CacheStrategySpec.DenseBoundaryCachingDetails) s.details;	
		Assert.a(details != null,"ERROR -- the Cache Strategy DenseNodeBoundary is too intricate to omit specifying details");
		
		int firstTableIndex = denseBoundaryTables.size();
		details.check();
		for (int j=0; j<details.nTables; j++) {
			float[] temp = new float[modelInfo.totalPositions+modelInfo.nSeqs];
			for (int k=0; k<modelInfo.totalPositions+modelInfo.nSeqs; k++) {
				temp[k] = 0;
			}
			denseBoundaryTables.add(temp);
		}
		
		ArrayList<DenseNodeBoundaryIndices> tempNodeBoundaries = new ArrayList<DenseNodeBoundaryIndices>();
		//for (int j=0; j<details.entries.size(); j++) {
		for (CacheStrategySpec.DenseBoundaryEntry be : details.entries) {
			DenseNodeBoundaryIndices db = new DenseNodeBoundaryIndices();
			
			int potential = be.potential; //details.potential[j];
			Assert.a(potential<modelInfo.nStates,"ERROR - Can't cache potential " + potential + " , which is a transition using the DENSE_NODE_BOUNDARY caching strategy.");
			//numFixedEvalIndices[potential]++;
			db.lookupTable = denseBoundaryTables.get(firstTableIndex + be.tableNum); //details.tableNum[j]);
			db.featureIndex = be.featureIndex; //details.featureIndex[j];
			db.potential = potential;
			db.rightPad = be.rightPad; //details.rightPad[j];
			db.leftPad = be.leftPad; //details.leftPad[j];
			denseNodeBoundaryEvalsList[potential].add(db);
			log.debug("Node boundary feature: "+db.featureIndex+" pot: "+db.potential+" table: "+db.lookupTable);
			tempNodeBoundaries.add(db);
		}
		Assert.a(tempNodeBoundaries.size() == details.entries.size());
		
		
		for (int seqNum=0; seqNum<data.size(); seqNum++) {
			TrainingSequence seq = data.get(seqNum);
			for (int pos=0; pos<seq.length(); pos++) {
				
				int cumulativePosition = modelInfo.seqOffsets[seqNum] + pos;
				
				for (int state=0; state< modelInfo.nStates; state++) {
					ArrayFeatureList result = new ArrayFeatureList(fm);				
					Assert.a(FeatureManagerNodeBoundaries.class.isInstance(fm1),"ERROR - to cache using DENSE NODE BOUNDARY, must be an instance of FeatureManagerNodeBoundaries, but this isn't.");
					FeatureManagerNode fmn = (FeatureManagerNode) fm1;
					fmn.evaluateNode(getSequence(seq, fmn),pos,state,result);
					
					if (result.isValid()) {	
						for (int j=0; j<result.size(); j++) {
							boolean found = false;
							for (DenseNodeBoundaryIndices db : tempNodeBoundaries ) {	
								if ((state==db.potential) && (result.getIndex(j)==db.featureIndex)) {
									int tx = modelInfo.seqOffsets[seqNum] + seqNum + pos;
									db.lookupTable[tx + 1] = db.lookupTable[tx] + (float) result.getValue(j);
									found = true;
									break;
								}
							}
							if(!found) {
								Assert.a(false, String.format("Feature Ix: %d State: %d not found. ", result.getIndex(j), state));
							}
						}
					} 
					
					if (!result.isValid()) {
						int ind = cumulativePosition*modelInfo.nPotentials + state;
						invalidTransitions[ind] = true;
					}
				}
			}
		}	
	}
	
	private void cacheLengthFeatureManagers() {
		// Iterate through each length for each state
		ArrayFeatureList featureList = new ArrayFeatureList(fm);
		for(CacheProcessor.StatePotentials stateInfo : modelInfo.statesWithLookback) {
			LengthOnlyEvaluation[] lengthForState = lengthTables[stateInfo.state];

			int maxLen = modelInfo.maxStateLengths[stateInfo.state];
			for(int i=0; i<maxLen; ++i) {
				LengthOnlyEvaluation eval = lengthForState[i];
				
				// Fill the evaluation object with the correct info for the state.
				featureList.clear();
				for(FeatureManagerNodeExplicitLength fm : (List<FeatureManagerNodeExplicitLength>) (List) lengthFMList) {
					fm.evaluateNodeLength(data.get(0),maxLen, i+1, stateInfo.state, featureList);
				}
				int size = featureList.currentSize;
				eval.featureIndex = new short[size]; 
				eval.value = new float[size];
				for(int j=0; j<size; ++j) {
					eval.featureIndex[j] = (short) featureList.indices[j];
					eval.value[j] = (float) featureList.values[j];
				}
			}
		}
	}
	
	public void checkConstraintsInTrainingData() {
		// Verify that there are no invalid transitions in the training data.
		List<TrainingSequence<?>> goodData = new ArrayList<TrainingSequence<?>>();

		for (int i=0; i<data.size(); ++i) {
			TrainingSequence seq = data.get(i);
			int seqOffset = modelInfo.seqOffsets[i];
			int lastState = -1;
			boolean validSequence = true;
			try {
				for (int pos=0; pos < seq.length(); pos++) {
					int trainingState = seq.getY(pos);
					int index = (seqOffset + pos) * modelInfo.nPotentials;
					Assert.a(!invalidTransitions[index + trainingState], "Seq: ",i," Pos: ", pos, " State: ",trainingState, " violates a constraint.");
					if(lastState != -1)
						Assert.a(!invalidTransitions[index + modelInfo.nStates + modelInfo.transitionIndex.getQuick(lastState, trainingState)], "Seq: ",i," Pos: ", pos, " Transition: ",lastState,"-",trainingState, " violates a constraint.");
					
					lastState = trainingState;
				}
			}
			catch(CheckException ex) {
				if(ignoreInvalidTrainingData) {
					log.warn("Discarding sequence "+i+" "+ex.getMessage());
					validSequence = false;
				}
				else {
					throw ex;
				}
			}
			if(validSequence) {
				goodData.add(seq);
			}
		}
		Assert.a(goodData.size() != 0, "All training sequences were invalid.");
		int diff = data.size() - goodData.size();
		if(diff != 0) {
			log.warn("Using "+goodData.size()+" training sequences.  Discarded "+diff+" because of constraint problems.");
			log.warn("Rebuilding the feature cache using the remaining good sequences.");
			// Redo the initialization with the new set of training data
			setTrainingData(fm, goodData);
			return;
		}
	}

	///////////////////////////////////////////////////////////////////////////////
	///// ABOVE: build the cache structure    BELOW: evaluate a position using the stored cache
	///////////////////////////////////////////////////////////////////////////////
	
	
	public void evaluatePosition(int seq, int pos) {
		// The job of this function is to update evals.index and evals.value based on all features and potentials at this location
		// Recall below the comments from another class about the information in evals:
		
		int cumulativePosition = modelInfo.seqOffsets[seq] + pos;		
		evaluateDenseFeatures(cumulativePosition);
		evaluateSparseFeatures(cumulativePosition); // only variable length
	}
	
	private void evaluateConstantFeatures() {
		Assert.a(!writtenYet, "Constant features can only be written once");
		for ( ConstantEvaluation ce : constEvals) {
			short fi    = ce.featureIndex;
			int pi      = ce.potential;
			int ci      = ce.evalIndex;
			float val   = ce.value;
			
			evals[pi].index[ci] = fi;
			evals[pi].value[ci] = val;			
		}
		writtenYet = true;
	}
	
	
	private void evaluateDenseFeatures(int cumulativePosition) {		
		for (DenseEvalIndices de : denseEvalsList) {			
			FeatureEvaluation fe  = de.evalTable;//   = densePotentialIndices[j];
			int               ei  = de.evalIndex;//denseEvalIndices[j];
			float[]           lut = de.lookupTable;//    = denseTableIndices[j];
			fe.index[ei] = de.featureIndex;
			fe.value[ei] = lut[cumulativePosition];
		}
	}
	
	
	private void evaluateSparseFeatures(int cumulativePosition) {
		for (int pot=0; pot<modelInfo.nPotentials; pot++) {
			currentEvalIndices[pot] = numFixedEvalIndices[pot];
		}
		int start = sparseStarts[cumulativePosition];
		int stop = sparseStarts[cumulativePosition+1];
		for (int j=start; j<stop; j++) {
			SparseEvaluation se = sparseEvals.get(j);
			short fi   = se.featureIndex;
			int pi     = se.potential;
			int ci     = currentEvalIndices[pi];
			currentEvalIndices[pi]++;
			float val  = se.value;
			
			evals[pi].index[ci] = fi;
			evals[pi].value[ci] = val;			
		}
		
		for (int pot=0; pot<modelInfo.nPotentials; pot++) {
			evals[pot].index[currentEvalIndices[pot]] = -1;
		}
	}
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////  ABOVE: Evaluate Markov features that don't depend on length   BELOW: evaluate semi-Markov features that depend on interval length
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	/**
	 *  The job of this function is to update the table "LengthFeatureEvaluation[][] lengthEvals"
	 *  
	 *  The first dimension is index by those states which have at least one explicit length node feature,
	 *  corresponding to the state of an interval which ends at position currently under consideration.
	 *  The size of this dimension is predetermined and is equal to modelInfo.statesWithLookback.length
	 *  
	 *  The second dimension is the length of the lookback, i.e. how many bases is
	 *  the length of interval currently being considered.  The size of this dimension is variable, and
	 *  after the last of which one inserts a LengthFeatureEvaluation whose lookback is -1
	 * 
	 *  A LengthFeatureEvaluation contains the lookback and a FeatureEvaluation for nodes, plus and
	 *  edgeEvaluation (which for now we set to null).  The FeatureEvaluations we have seen before;
	 *  they are comprised of index and value arrays, which are of variable length,
	 *  and after last entry you put a -1 in the index array.
	 * 
	 */
	
	public void evaluateSegmentsEndingAt(int seq, int pos) {
		
		int seqOffset = modelInfo.seqOffsets[seq];
		int overallPosition = seqOffset+pos;
		int tx1 = modelInfo.seqOffsets[seq]+pos + seq + 1;  // ending position on subtraction lookup table, remember each sequence needs one extra position of padding.
		
		CacheProcessor.StatePotentials[] statesWithLookback = modelInfo.statesWithLookback;
		int nSemiMarkovStates = statesWithLookback.length;

		int seqLen = data.get(seq).length();
		int invalidIndex = overallPosition*modelInfo.nPotentials; 
		for (int stateIx=0; stateIx < nSemiMarkovStates; stateIx++) {
			CacheProcessor.StatePotentials statePotentials = statesWithLookback[stateIx];
			LengthFeatureEvaluation[] lookbackEvals = lengthEvals[stateIx];

			if(!checkExit(invalidIndex, pos, seqLen, statePotentials.state)) {
				lookbackEvals[0].lookback = -1;
				continue;
			}
			
			// For say an interval of exon1 ending at specified position, what lookback need we provide?
			
			// Well, lets say for instance that exons have a minimum length of 50 and maximum length of 200.
			// Then I want to start at current pposition (lookback=0) and walk backwards, breaking if I ever
			// encounter an invalidated node or edge, or when I reach the maximum allowable length for exon.
			// If not, then as soon as I get to lookback=50, I
			// start taking notes for every valid beginning (a valid edge leading into beginning of current interval)
			
			// Now for one of these, I look back, and for every valid transition into I record the lookback, the
			// relevant interval node features, both index and value, by subtracting some numbers at an offset.
			
			int nLookbacks = 0;
			int minLength = minStateLengths[statesWithLookback[stateIx].state];
			int maxLength = maxStateLengths[statesWithLookback[stateIx].state];

			for (int lookback = 0; lookback < maxLength; lookback++) {
				int firstPosIndex = (overallPosition-lookback)*modelInfo.nPotentials;
				if (lookback > pos || invalidTransitions[firstPosIndex + statePotentials.state]) { break; }
				if (lookback+1 < minLength) { continue; }
				
				int prevPos = pos - lookback - 1;
				boolean validEntry = false;
				if(prevPos == -1) {
					validEntry = true;
				}
				else {
					for (int pot : statePotentials.potentials) {
						int entryIndex = (seqOffset+prevPos+1)*modelInfo.nPotentials;
						if (modelInfo.selfTransitions[statePotentials.state]+modelInfo.nStates != pot && !invalidTransitions[entryIndex + pot]) {
							validEntry = true;
							break; 
						}
					}
				}
				
				if(validEntry) {
					int nEvals = 0;
					LengthFeatureEvaluation lengthEval = lookbackEvals[nLookbacks];
					lengthEval.lookback  = (short) lookback;
					lengthEval.edgeEvals = null;
					
					FeatureEvaluation nodeEval = lengthEval.nodeEval;

						// for this potential, there is a lookup table.  It includes the offsets you need for subtraction etc.
					for (DenseNodeBoundaryIndices db : denseNodeBoundaryEvalsList[statePotentials.state]) {
						int index    = db.featureIndex;
						float[] lut  = db.lookupTable;
						int rightPad = db.rightPad;
						int leftPad  = db.leftPad;
						
						float val = lut[tx1 - rightPad ] - lut[tx1-lookback-1 + leftPad];
						
						nodeEval.index[nEvals] = (short) index;
						nodeEval.value[nEvals] = val;
						nEvals++;
					}

					// Add in length evaluations
					LengthOnlyEvaluation lengthOnlyEval = lengthTables[statePotentials.state][lookback];
					int size = lengthOnlyEval.featureIndex.length;
					System.arraycopy(lengthOnlyEval.featureIndex, 0, nodeEval.index, nEvals, size);
					System.arraycopy(lengthOnlyEval.value, 0, nodeEval.value, nEvals, size);
					nEvals += size;
					
					nodeEval.index[nEvals] = -1;						
					nLookbacks++;
				}
			}
			lookbackEvals[nLookbacks].lookback = -1;
		}
	}
	
	/** Checks if there is a valid transition out of a node */
	boolean checkExit(int positionIndex, int pos, int seqLen, int state) {
		if(invalidTransitions[positionIndex+state])
			return false;

		// This requires that we be in the last position or that there is a valid transition out.
		if(pos == seqLen-1) {
			return true;
		}
		
		boolean wayOut = false;
		int nextPosIndex = positionIndex +  modelInfo.nPotentials;
		for(byte pot : modelInfo.exitTransitions[state]) {
			if(modelInfo.selfTransitions[state]+modelInfo.nStates != pot && !invalidTransitions[nextPosIndex + pot]) {
				wayOut = true;
				break;
			}
		}
		return wayOut; 
	}
	
	/////////////////////////////////////////////////////////////////////////////////////////////////
	
	
	public boolean[] getInvalidTransitions() {
		return invalidTransitions;
	}
	
	public void setSemiMarkovSetup(SemiMarkovSetup setup) {
		maxStateLengths = setup.getMaxLengths();
		minStateLengths = setup.getMinLengths();
		ignoreSemiMarkovSelfTransitions = setup.isIgnoreSemiMarkovSelfTransitions();
	}
		
	public boolean isAllPaths() {
		return allPaths;
	}
	
	
	public void setAllPaths(boolean allPaths) {
		this.allPaths = allPaths;
	}

	/**
	 * @return Returns the lookbackArrayFeatureSize.
	 */
	public int getLookbackArrayFeatureSize() {
		return lookbackArrayFeatureSize;
	}

	/**
	 * @param lookbackArrayFeatureSize The lookbackArrayFeatureSize to set.
	 */
	public void setLookbackArrayFeatureSize(int lookbackArrayFeatureSize) {
		this.lookbackArrayFeatureSize = lookbackArrayFeatureSize;
	}

	/**
	 * @return Returns the lookbackArraySize.
	 */
	public int getLookbackArraySize() {
		return lookbackArraySize;
	}

	/**
	 * @param lookbackArraySize The lookbackArraySize to set.
	 */
	public void setLookbackArraySize(int lookbackArraySize) {
		this.lookbackArraySize = lookbackArraySize;
	}

	public boolean isIgnoreInvalidTrainingData() {
		return ignoreInvalidTrainingData;
	}

	public void setIgnoreInvalidTrainingData(boolean ignoreInvalidTrainingData) {
		this.ignoreInvalidTrainingData = ignoreInvalidTrainingData;
	}
	
	private final static void debug(Object msg){
		if(LOGGER.isDebugEnabled()){
			LOGGER.debug(msg.toString());
		}
	}
	
	private final static void error(Object msg){
		
			LOGGER.error(msg.toString());
		
	}
}
