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

public class PfamGenic extends AbstractFeatureManager<CompositeInput> implements FeatureManagerNode<CompositeInput> {
	private static final long serialVersionUID = -7659288739348604129L;
	private static final Log log = LogFactory.getLog(PfamGenic.class);
	boolean debug = log.isDebugEnabled();
	
	/* Contains 1 features:
	 *    f returns 1 if either of two conditions below and 0 otherwise:
	 *      a) y_i=intron1,intron2,intron3 and pest(i+1) = 2 [intron only]
	 *      b) y_i= intron1m,intron2m,intron3m and mest(i+1) = 2 [intron only]
	 */
	
	int startIx;  // The index of the first feature managed by this FeatureManager
	ModelManager model;
	
	boolean[] plusGenic;
	boolean[] minusGenic;
	
	
	public PfamGenic() {
	}

	public int getNumFeatures() {
		return 1;
	}	
	
	public String getFeatureName(int featureIndex) {
		return "PfamGenic";
	}

	
	public void evaluateNode(InputSequence<? extends CompositeInput> seq, int pos, int state, FeatureList result) {
		if(pos == seq.length()-1) {
			return;
		}		

		InputSequence<Integer>  ppfam = (InputSequence<Integer>) seq.getComponent("ppfam");		
		InputSequence<Integer>  mpfam = (InputSequence<Integer>) seq.getComponent("mpfam");
		
		boolean plusPfam = (ppfam.getX(pos) > 0);
		boolean minusPfam = (mpfam.getX(pos) > 0);		

		if (plusGenic[state] && plusPfam) { result.addFeature(startIx, 1); }
		if (minusGenic[state] && minusPfam) { result.addFeature(startIx, 1); }		
	}


	public void train(int startingIndex, ModelManager modelInfo, List<? extends TrainingSequence<? extends CompositeInput>> data) {
		startIx = startingIndex;
		model = modelInfo;

		int nStates = model.getNumStates();
		
		plusGenic = new boolean[nStates];
		for (int j=0; j<nStates; j++) { plusGenic[j] = false; }
		plusGenic[model.getStateIndex("exon1")] = true;
		plusGenic[model.getStateIndex("exon2")] = true;
		plusGenic[model.getStateIndex("exon3")] = true;
		plusGenic[model.getStateIndex("intron1")] = true;
		plusGenic[model.getStateIndex("intron2")] = true;
		plusGenic[model.getStateIndex("intron3")] = true;
		
		minusGenic = new boolean[nStates]; 
		for (int j=0; j<nStates; j++) { minusGenic[j] = false; }		
		minusGenic[model.getStateIndex("exon1m")] = true;
		minusGenic[model.getStateIndex("exon2m")] = true;
		minusGenic[model.getStateIndex("exon3m")] = true;
		minusGenic[model.getStateIndex("intron1m")] = true;
		minusGenic[model.getStateIndex("intron2m")] = true;
		minusGenic[model.getStateIndex("intron3m")] = true;		
	}
	@Override
	public CacheStrategySpec getCacheStrategy() {
		return new CacheStrategySpec(CacheStrategy.UNSPECIFIED);
	}
}

