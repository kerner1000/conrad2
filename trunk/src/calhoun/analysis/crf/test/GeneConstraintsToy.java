package calhoun.analysis.crf.test;

import java.util.List;

import calhoun.analysis.crf.AbstractFeatureManager;
import calhoun.analysis.crf.FeatureList;
import calhoun.analysis.crf.FeatureManagerEdge;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.analysis.crf.io.TrainingSequence;

/** Implements basic constraints on gene calls.
 * 
 * 1) Intergenic - start must occur at ATG
 * 2) Splice sites must be canonical GT/AG or GC/AG
 * 3) Exon-stop must be followed by a start codon
 */

// NOTE: I don't think these constraints look at frame, so I recommend using a version of this class adapted to the specific model yu're considering.  JPV 20060629

public class GeneConstraintsToy extends AbstractFeatureManager<Character> implements FeatureManagerEdge<Character> {
	
	private static final long serialVersionUID = -3753476830756229273L;

	public String getFeatureName(int featureIndex) {
		return "Gene constraints toy";
	}

	/** This is a constraint class, so we don't return features */
	public int getNumFeatures() {
		return 0;
	}

	/** Set up the matrix
	 * Depends on states starting with the words 'intergenic, intron, and exon'.  Also depends on the negative strand states ending in m.
	 */
	public void train(int startingIndex, ModelManager modelInfo, List<? extends TrainingSequence<? extends Character>> data) {

	}
	
	public void evaluateEdge(InputSequence<? extends Character> seq, int pos, int prevState, int state, FeatureList result) {
		if (pos==500) {
			if  ( (state != 1) || (prevState!=0) ) {
				//System.out.println("Just triggered the constraint in GeneConstraintsToy");
				result.invalidate();
			}
		}
	}

}
