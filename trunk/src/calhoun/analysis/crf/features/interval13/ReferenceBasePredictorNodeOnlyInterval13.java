package calhoun.analysis.crf.features.interval13;

import calhoun.analysis.crf.CacheStrategySpec;
import calhoun.analysis.crf.FeatureManagerNode;
import calhoun.analysis.crf.CacheStrategySpec.CacheStrategy;


public class ReferenceBasePredictorNodeOnlyInterval13 extends ReferenceBasePredictorInterval13Base implements FeatureManagerNode<Character> {
	
	private static final long serialVersionUID = 6351525014233239165L;

	public ReferenceBasePredictorNodeOnlyInterval13() {
	}

	@Override
	public CacheStrategySpec getCacheStrategy() {

		CacheStrategySpec css = new CacheStrategySpec(CacheStrategy.DENSE);		

		return css;
	}

}
