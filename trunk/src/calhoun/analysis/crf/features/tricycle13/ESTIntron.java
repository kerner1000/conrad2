package calhoun.analysis.crf.features.tricycle13;

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

public class ESTIntron extends AbstractFeatureManager<CompositeInput> implements FeatureManagerNode<CompositeInput> {
	private static final long serialVersionUID = -7659288739348604129L;
	private static final Log log = LogFactory.getLog(ESTIntron.class);
	boolean debug = log.isDebugEnabled();
	
	/* Contains 1 features:
	 *    f returns 1 if either of two conditions below and 0 otherwise:
	 *      a) y_i=intron1,intron2,intron3 and pest(i+1) = 2 [intron only]
	 *      b) y_i= intron1m,intron2m,intron3m and mest(i+1) = 2 [intron only]
	 */
	
	int startIx;  // The index of the first feature managed by this FeatureManager
	ModelManager model;
	
	boolean[] plusIntronState;
	boolean[] minusIntronState;
	
	
	public ESTIntron() {
	}

	public int getNumFeatures() {
		return 1;
	}	
	
	public String getFeatureName(int featureIndex) {
		return "ESTIntron";
	}
	
	public void evaluateNode(InputSequence<? extends CompositeInput> seq, int pos, int state, FeatureList result) {
		if(pos == seq.length()-1) {
			return;
		}		

		InputSequence<Integer>  pest = (InputSequence<Integer>) seq.getComponent("pest");		
		InputSequence<Integer>  mest = (InputSequence<Integer>) seq.getComponent("mest");
		
		boolean plusEstIntron = (pest.getX(pos+1) == 2);
		boolean minusEstIntron = (mest.getX(pos+1) == 2);		

		if (plusIntronState[state] && plusEstIntron) { result.addFeature(startIx, 1); }
		if (minusIntronState[state] && minusEstIntron) { result.addFeature(startIx, 1); }		
	}


	public void train(int startingIndex, ModelManager modelInfo, List<? extends TrainingSequence<? extends CompositeInput>> data) {
		startIx = startingIndex;
		model = modelInfo;

		int nStates = model.getNumStates();
		
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
	@Override
	public CacheStrategySpec getCacheStrategy() {
		return new CacheStrategySpec(CacheStrategy.UNSPECIFIED);
	}
}

