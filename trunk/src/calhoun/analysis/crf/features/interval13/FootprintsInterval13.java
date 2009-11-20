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
import calhoun.analysis.crf.io.InputSequence;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.analysis.crf.io.MultipleAlignmentInputSequence.MultipleAlignmentColumn;
import calhoun.seq.KmerHasher;
import calhoun.util.Assert;

public class FootprintsInterval13 extends AbstractFeatureManager<MultipleAlignmentColumn> implements FeatureManagerNode<MultipleAlignmentColumn> {
	private static final long serialVersionUID = -7659288739348604129L;
	private static final Log log = LogFactory.getLog(FootprintsInterval13.class);
	boolean debug = log.isDebugEnabled();
	
	/* Features are the conjunction of "species X is present in multiple alignment" with hidden state is "exonic, intronic, intergenic"
	 * Is the number of features allowed to depend on the number of species inmultiple alignment??
	 */

	List<String> speciesNames;
	int startIx;  
	ModelManager model;
	KmerHasher h = new KmerHasher(KmerHasher.ACGTN,1);
	
	int maxSeqLength;
	
	int nFeatures = -1;
	
	Boolean[] isStateCoding, isStateIntronic, isStateIntergenic;

	
	public FootprintsInterval13() {	
	}

	public int getNumFeatures() {
		return nFeatures;
	}	
	
	public String getFeatureName(int featureIndex) {
		String[] type = new String[] { "intergenic", "exonic", "intronic"};
		int raw = featureIndex - startIx;
		Assert.a(raw<nFeatures);
		if(speciesNames == null) {
			return "Species "+((raw/3) + 1) + " "+type[raw%3]+" footprint";
		}
		return speciesNames.get((raw/3) + 1) + " "+type[raw%3]+" footprint";
	}
	
	
	public void evaluateNode(InputSequence<? extends MultipleAlignmentColumn> seq, int pos, int state, FeatureList result) {
		MultipleAlignmentColumn mac = seq.getX(pos);
		for (int species = 1; species<mac.numSpecies(); species++) {
			if (mac.nucleotide(species) == '-') continue;
			
			if (isStateIntergenic[state]) { result.addFeature(startIx+((species-1)*3+0), 1.0); }
			if (isStateCoding[state])     { result.addFeature(startIx+((species-1)*3+1), 1.0); }
			if (isStateIntronic[state])   { result.addFeature(startIx+((species-1)*3+2), 1.0); }
		}
	}


	public void train(int startingIndex, ModelManager modelInfo, List<? extends TrainingSequence<? extends MultipleAlignmentColumn>> data) {
		TrainingSequence<? extends MultipleAlignmentColumn> seq = data.get(0);
		speciesNames = seq.getX(0).getMultipleAlignment().getSpeciesNames();
		
		startIx = startingIndex;
		model = modelInfo;
		int nStates = model.getNumStates();

		nFeatures = 3*(data.get(0).getX(0).numSpecies()-1);  // Assumes this is the same for all alignments
						
		isStateCoding = new Boolean[nStates];       for (int j=0; j<nStates; j++) { isStateCoding[j] = false; }
		isStateCoding[1] = true;
		isStateCoding[2] = true;
		isStateCoding[3] = true;
		isStateCoding[7] = true;
		isStateCoding[8] = true;
		isStateCoding[9] = true;		

		isStateIntronic = new Boolean[nStates];     for (int j=0; j<nStates; j++) { isStateIntronic[j] = false; }
		isStateIntronic[4] = true;
		isStateIntronic[5] = true;
		isStateIntronic[6] = true;
		isStateIntronic[10] = true;
		isStateIntronic[11] = true;
		isStateIntronic[12] = true;

		isStateIntergenic = new Boolean[nStates];   for (int j=0; j<nStates; j++) { isStateIntergenic[j] = false; }
		isStateIntergenic[0] = true;
		
	}
	@Override
	public CacheStrategySpec getCacheStrategy() {
		return new CacheStrategySpec(CacheStrategy.DENSE);
	}
	
}

