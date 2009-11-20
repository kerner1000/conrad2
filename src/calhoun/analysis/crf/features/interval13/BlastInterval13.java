package calhoun.analysis.crf.features.interval13;

import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.AbstractFeatureManager;
import calhoun.analysis.crf.CacheStrategySpec;
import calhoun.analysis.crf.FeatureList;
import calhoun.analysis.crf.FeatureManagerNode;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.CacheStrategySpec.CacheStrategy;
import calhoun.analysis.crf.io.CompositeInput;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.analysis.crf.io.TrainingSequence;

public class BlastInterval13 extends AbstractFeatureManager<CompositeInput> implements FeatureManagerNode<CompositeInput> {
	private static final long serialVersionUID = -7659288739348604129L;
	private static final Log log = LogFactory.getLog(BlastInterval13.class);
	boolean debug = log.isDebugEnabled();
	
	int startIx;  // The index of the first feature managed by this FeatureManager
	ModelManager model;
	
	boolean[] intergenicState;
	boolean[] plusExonState;
	boolean[] minusExonState;
	boolean[] plusIntronState;
	boolean[] minusIntronState;	
	
	public BlastInterval13() {
	}

	public int getNumFeatures() {
		return 3;
	}	
	
	public String getFeatureName(int featureIndex) {
		String[] types = new String[] {"exon", "intron", "Intergenic","a"};
		return "Blast "+types[featureIndex-startIx];
	}

	public void evaluateNode(InputSequence<? extends CompositeInput> seq, int pos, int state, FeatureList result) {
		if(pos == seq.length()-1) {
			return;
		}		

		InputSequence<Integer>  pblast = (InputSequence<Integer>) seq.getComponent("pblast");
		InputSequence<Integer>  mblast = (InputSequence<Integer>) seq.getComponent("mblast");

		int plusEst  = pblast.getX(pos+1);
		int minusEst = mblast.getX(pos+1);
		// 0 - no data
		// 1 - blast cluster hit here on this strand

		if ((plusExonState[state] && (plusEst==1)) || (minusExonState[state] && (minusEst==1))) { 
			result.addFeature(startIx, 1); 
		}
		
		if ((plusIntronState[state] && (plusEst==1)) || (minusIntronState[state] && (minusEst==1))) { 
			result.addFeature(startIx+1, 1); 
		}

		if (intergenicState[state] && (plusEst==1 || minusEst==1)) { 
			result.addFeature(startIx+2, 1); 
		}
}
		
	public void train(int startingIndex, ModelManager modelInfo, List<? extends TrainingSequence<? extends CompositeInput>> data) {
		startIx = startingIndex;
		model = modelInfo;

		int nStates = model.getNumStates();
		
		plusExonState = new boolean[nStates];
		plusExonState[model.getStateIndex("exon0")] = true;
		plusExonState[model.getStateIndex("exon1")] = true;
		plusExonState[model.getStateIndex("exon2")] = true;
		
		minusExonState = new boolean[nStates]; 
		minusExonState[model.getStateIndex("exon0m")] = true;
		minusExonState[model.getStateIndex("exon1m")] = true;
		minusExonState[model.getStateIndex("exon2m")] = true;		
		
		plusIntronState = new boolean[nStates];
		plusIntronState[model.getStateIndex("intron0")] = true;
		plusIntronState[model.getStateIndex("intron1")] = true;
		plusIntronState[model.getStateIndex("intron2")] = true;
		
		minusIntronState = new boolean[nStates]; 
		minusIntronState[model.getStateIndex("intron0m")] = true;
		minusIntronState[model.getStateIndex("intron1m")] = true;
		minusIntronState[model.getStateIndex("intron2m")] = true;	

		intergenicState = new boolean[nStates];
		intergenicState[model.getStateIndex("intergenic")] = true;
	}

	@Override
	public CacheStrategySpec getCacheStrategy() {
		return new CacheStrategySpec(CacheStrategy.DENSE);
	}
}

