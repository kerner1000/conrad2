package calhoun.analysis.crf.features.interval13;

import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.AbstractFeatureManager;
import calhoun.analysis.crf.CacheStrategySpec;
import calhoun.analysis.crf.FeatureList;
import calhoun.analysis.crf.FeatureManagerEdge;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.CacheStrategySpec.CacheStrategy;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.util.Assert;


public class StateTransitionsInterval13 extends AbstractFeatureManager<Character> implements FeatureManagerEdge<Character> {
	private static final long serialVersionUID = -7745853849576425025L;
	private static final Log log = LogFactory.getLog(StateTransitionsInterval13.class);

	private int startIx;
	private double intronProb;
	private double endProb;
	
	public String getFeatureName(int featureIndex) {
		return "State transitions";
	}

	public int getNumFeatures() {
		return 1;
	}

	public void train(int startingIndex, ModelManager modelInfo, List<? extends TrainingSequence<? extends Character>> data) {
		Interval13Tools.verify(modelInfo);
		startIx = startingIndex;
		
		// Get the average # of exons per gene from training data
		int intronCount = 0;
		int geneCount = 0;
		for(TrainingSequence<?> seq : data) {
			int[] y = seq.getY();
			int prevState = y[0];
			for(int i=1; i<y.length; ++i) {
				int state = y[i];
				switch(Interval13Tools.edgeConstraints[prevState*Interval13Tools.numStates + state]) {
					case PDON:
					case MACC:
						++intronCount;
						break;
					case PSTOP:
					case MSTART:
						++geneCount;
						break;
					default:
				}
				prevState = state;
			}
		}
		
		double avgExonCount = (intronCount+geneCount)/((float)geneCount);
		endProb = Math.log(1/avgExonCount);
		intronProb = Math.log(1 - 1/avgExonCount);
		log.warn(String.format("%d genes, %d introns, %.2f exons/gene", geneCount, intronCount, avgExonCount));
	}
	
	public void evaluateEdge(InputSequence<? extends Character> seq, int pos, int prevState, int state, FeatureList result) {
	
		// There's really only one parameter below: the average number of exons in a gene.
		
		switch(Interval13Tools.edgeConstraints[prevState*Interval13Tools.numStates + state]) {
			case NONE:
			case PACC:
			case MDON:
			case PCODE: // redundant with node invalidation below
			case MCODE: // redundant iwth node evaluation below
				break;
			case PSTART:
			case MSTOP:
				result.addFeature(startIx,Math.log(0.5));
				break;
			case PDON:
			case MACC:
				result.addFeature(startIx, intronProb);
				break;
			case PSTOP:
			case MSTART:
				result.addFeature(startIx, endProb);
				break;
			case NEVER:
			default:
				Assert.a(false);
		}
	}

	
	@Override
	public CacheStrategySpec getCacheStrategy() {
		return new CacheStrategySpec(CacheStrategy.CONSTANT);
	}

	/**
	 * @return Returns the endProb.
	 */
	public double getEndProb() {
		return endProb;
	}

	/**
	 * @return Returns the intronProb.
	 */
	public double getIntronProb() {
		return intronProb;
	}

}
