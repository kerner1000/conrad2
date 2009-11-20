package calhoun.analysis.crf.features.tricycle13;

import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.AbstractFeatureManager;
import calhoun.analysis.crf.FeatureList;
import calhoun.analysis.crf.FeatureManagerNode;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.io.CompositeInput;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.analysis.crf.io.TrainingSequence;

public class ESTExon extends AbstractFeatureManager<CompositeInput> implements FeatureManagerNode<CompositeInput> {
	private static final long serialVersionUID = -7659288739348604129L;
	private static final Log log = LogFactory.getLog(ESTExon.class);
	boolean debug = log.isDebugEnabled();
	

	
	int startIx;  // The index of the first feature managed by this FeatureManager
	ModelManager model;
	
	boolean[] plusExonState;
	boolean[] minusExonState;
	boolean[] plusIntronState;
	boolean[] minusIntronState;	
	
	public ESTExon() {
	}

	public int getNumFeatures() {
		return 4;
	}	
	
	public String getFeatureName(int featureIndex) {
		return "ESTExon";
	}
	
	public void evaluateNode(InputSequence<? extends CompositeInput> seq, int pos, int state, FeatureList result) {
		if(pos == seq.length()-1) {
			return;
		}		

		InputSequence<Integer>  pest = (InputSequence<Integer>) seq.getComponent("pest");		
		InputSequence<Integer>  mest = (InputSequence<Integer>) seq.getComponent("mest");
		
		int plusEst  = pest.getX(pos+1);
		int minusEst = mest.getX(pos+1);
		// 0 - no data
		// 1 - exon only
		// 2 - intron only
		// 3 - mixed

		if (plusExonState[state] && (plusEst==1)) { result.addFeature(startIx, 1); }
		if (minusExonState[state] && (minusEst==1)) { result.addFeature(startIx, 1); }

		if (plusExonState[state] && (plusEst==3)) { result.addFeature(startIx+1, 1); }
		if (minusExonState[state] && (minusEst==3)) { result.addFeature(startIx+1, 1); }
		
		if (plusIntronState[state] && (plusEst==2)) { result.addFeature(startIx+2, 1); }
		if (minusIntronState[state] && (minusEst==2)) { result.addFeature(startIx+2, 1); }

		if (plusIntronState[state] && (plusEst==3)) { result.addFeature(startIx+3, 1); }
		if (minusIntronState[state] && (minusEst==3)) { result.addFeature(startIx+3, 1); }		
		
	}


	public void train(int startingIndex, ModelManager modelInfo, List<? extends TrainingSequence<? extends CompositeInput>> data) {
		startIx = startingIndex;
		model = modelInfo;

		int nStates = model.getNumStates();
		
		plusExonState = new boolean[nStates];
		for (int j=0; j<nStates; j++) { plusExonState[j] = false; }
		plusExonState[model.getStateIndex("exon1")] = true;
		plusExonState[model.getStateIndex("exon2")] = true;
		plusExonState[model.getStateIndex("exon3")] = true;
		
		minusExonState = new boolean[nStates]; 
		for (int j=0; j<nStates; j++) { minusExonState[j] = false; }		
		minusExonState[model.getStateIndex("exon1m")] = true;
		minusExonState[model.getStateIndex("exon2m")] = true;
		minusExonState[model.getStateIndex("exon3m")] = true;		
		
		plusIntronState = new boolean[nStates];
		for (int j=0; j<nStates; j++) { plusIntronState[j] = false; }
		plusIntronState[model.getStateIndex("intron1")] = true;
		plusIntronState[model.getStateIndex("intron2")] = true;
		plusIntronState[model.getStateIndex("intron3")] = true;
		
		minusIntronState = new boolean[nStates]; 
		for (int j=0; j<nStates; j++) { minusIntronState[j] = false; }		
		minusIntronState[model.getStateIndex("intron1m")] = true;
		minusIntronState[model.getStateIndex("intron2m")] = true;
		minusIntronState[model.getStateIndex("intron3m")] = true;	
	}
}

