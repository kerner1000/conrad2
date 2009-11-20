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

public class PfamPhase extends AbstractFeatureManager<CompositeInput> implements FeatureManagerNode<CompositeInput> {
	private static final long serialVersionUID = -7659288739348604129L;
	private static final Log log = LogFactory.getLog(PfamPhase.class);
	boolean debug = log.isDebugEnabled();
	
	/* Contains 1 features:
	 *    f returns 1 if either of two conditions below and 0 otherwise:
	 *      a) y_i=intron1,intron2,intron3 and pest(i+1) = 2 [intron only]
	 *      b) y_i= intron1m,intron2m,intron3m and mest(i+1) = 2 [intron only]
	 */
	
	int startIx;  // The index of the first feature managed by this FeatureManager
	ModelManager model;
	
	int exon1S;
	int exon2S;
	int exon3S;
	int exon1mS;
	int exon2mS;
	int exon3mS;
	
	public PfamPhase() {
	}

	public int getNumFeatures() {
		return 1;
	}	
	
	public String getFeatureName(int featureIndex) {
		return "PfamPhase";
	}

	
	public void evaluateNode(InputSequence<? extends CompositeInput> seq, int pos, int state, FeatureList result) {
		if(pos == seq.length()-1) {
			return;
		}		

		InputSequence<Integer>  ppfam = (InputSequence<Integer>) seq.getComponent("ppfam");		
		InputSequence<Integer>  mpfam = (InputSequence<Integer>) seq.getComponent("mpfam");
		
		if ((state==exon1S) && (ppfam.getX(pos)==1))  { result.addFeature(startIx, 1); }
		if ((state==exon2S) && (ppfam.getX(pos)==2))  { result.addFeature(startIx, 1); }
		if ((state==exon3S) && (ppfam.getX(pos)==3))  { result.addFeature(startIx, 1); }

		if ((state==exon1mS) && (mpfam.getX(pos)==1))  { result.addFeature(startIx, 1); }
		if ((state==exon2mS) && (mpfam.getX(pos)==2))  { result.addFeature(startIx, 1); }
		if ((state==exon3mS) && (mpfam.getX(pos)==3))  { result.addFeature(startIx, 1); }
	}


	public void train(int startingIndex, ModelManager modelInfo, List<? extends TrainingSequence<? extends CompositeInput>> data) {
		startIx = startingIndex;
		model = modelInfo;

		exon1S = model.getStateIndex("exon1");
		exon2S = model.getStateIndex("exon2");
		exon3S = model.getStateIndex("exon3");

		exon1mS = model.getStateIndex("exon1m");
		exon2mS = model.getStateIndex("exon2m");
		exon3mS = model.getStateIndex("exon3m");
	}
	@Override
	public CacheStrategySpec getCacheStrategy() {
		return new CacheStrategySpec(CacheStrategy.UNSPECIFIED);
	}
}

