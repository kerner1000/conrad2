package calhoun.analysis.crf.features.interval29;

import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.AbstractFeatureManager;
import calhoun.analysis.crf.CacheStrategySpec;
import calhoun.analysis.crf.FeatureList;
import calhoun.analysis.crf.FeatureManagerEdge;
import calhoun.analysis.crf.FeatureManagerNode;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.CacheStrategySpec.CacheStrategy;
import calhoun.analysis.crf.io.CompositeInput;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.analysis.crf.io.TrainingSequence;

public class ESTInterval29 extends AbstractFeatureManager<CompositeInput> implements FeatureManagerNode<CompositeInput> , FeatureManagerEdge<CompositeInput> {
	private static final long serialVersionUID = -7659288739348604129L;
	private static final Log log = LogFactory.getLog(ESTInterval29.class);
	boolean debug = log.isDebugEnabled();
	

	
	int startIx;  // The index of the first feature managed by this FeatureManager
	ModelManager model;
	
	boolean[] plusExonState;
	boolean[] minusExonState;
	boolean[] plusIntronState;
	boolean[] minusIntronState;	
	
	public ESTInterval29() {
	}

	public int getNumFeatures() {
		return 6;
	}	
	
	public String getFeatureName(int featureIndex) {
		return "ESTEdge";
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

	public void evaluateEdge(InputSequence<? extends CompositeInput> seq, int pos, int prevState, int state, FeatureList result) {
		if( (pos == seq.length()-1) || (pos == 0) ) {
			return;
		}		
		
		InputSequence<Integer>  pest = (InputSequence<Integer>) seq.getComponent("pest");
		InputSequence<Integer>  mest = (InputSequence<Integer>) seq.getComponent("mest");

		int plusE  = pest.getX(pos+1); int plusEp = pest.getX(pos);
		int minusE = mest.getX(pos+1); int minusEp = mest.getX(pos);
		
		Boolean plusacc = ((plusEp==2) && (plusE==1));
		Boolean plusdon = ((plusEp==1) && (plusE==2));
		Boolean minusacc = ((minusEp==1) && (minusE==2));
		Boolean minusdon = ((minusEp==2) && (minusE==1));
		
		if (plusacc  && plusExonState[state] && plusIntronState[prevState]) { result.addFeature(startIx+4,1); }
		if (minusacc && minusExonState[prevState] && minusIntronState[state]) { result.addFeature(startIx+4,1); }
		
		if (plusdon && plusIntronState[state] && plusExonState[prevState]) { result.addFeature(startIx+5,1); }
		if (minusdon && plusIntronState[prevState] && plusExonState[state]) { result.addFeature(startIx+5,1); }

	}
		
	public void train(int startingIndex, ModelManager modelInfo, List<? extends TrainingSequence<? extends CompositeInput>> data) {
		startIx = startingIndex;
		model = modelInfo;

		int nStates = model.getNumStates();
		
		plusExonState = new boolean[nStates];
		for (int j=0; j<nStates; j++) { plusExonState[j] = false; }
		plusExonState[model.getStateIndex("exon0")] = true;
		plusExonState[model.getStateIndex("exon1")] = true;
		plusExonState[model.getStateIndex("exon2")] = true;
		
		minusExonState = new boolean[nStates]; 
		for (int j=0; j<nStates; j++) { minusExonState[j] = false; }		
		minusExonState[model.getStateIndex("exon0m")] = true;
		minusExonState[model.getStateIndex("exon1m")] = true;
		minusExonState[model.getStateIndex("exon2m")] = true;		
		
		plusIntronState = new boolean[nStates];
		for (int j=0; j<nStates; j++) { plusIntronState[j] = false; }
		plusIntronState[model.getStateIndex("e-i0")] = true;
		plusIntronState[model.getStateIndex("e-i1")] = true;
		plusIntronState[model.getStateIndex("e-i2")] = true;
		plusIntronState[model.getStateIndex("intron0")] = true;
		plusIntronState[model.getStateIndex("intron1")] = true;
		plusIntronState[model.getStateIndex("intron2")] = true;
		plusIntronState[model.getStateIndex("i-e0")] = true;
		plusIntronState[model.getStateIndex("i-e1")] = true;
		plusIntronState[model.getStateIndex("i-e2")] = true;		
		
		minusIntronState = new boolean[nStates]; 
		for (int j=0; j<nStates; j++) { minusIntronState[j] = false; }		
		minusIntronState[model.getStateIndex("em-i0m")] = true;
		minusIntronState[model.getStateIndex("em-i1m")] = true;
		minusIntronState[model.getStateIndex("em-i2m")] = true;	
		minusIntronState[model.getStateIndex("intron0m")] = true;
		minusIntronState[model.getStateIndex("intron1m")] = true;
		minusIntronState[model.getStateIndex("intron2m")] = true;
		minusIntronState[model.getStateIndex("im-e0m")] = true;
		minusIntronState[model.getStateIndex("im-e1m")] = true;
		minusIntronState[model.getStateIndex("im-e2m")] = true;			
	}

	public CacheStrategySpec getCacheStrategy() {
		return new CacheStrategySpec(CacheStrategy.DENSE);
	}

}

