package calhoun.analysis.crf.features.interval29;

import calhoun.analysis.crf.CacheStrategySpec;
import calhoun.analysis.crf.FeatureManagerNode;
import calhoun.analysis.crf.CacheStrategySpec.CacheStrategy;


public class ReferenceBasePredictorNodeOnlyInterval29 extends ReferenceBasePredictorInterval29Base implements FeatureManagerNode<Character> {
	private static final long serialVersionUID = 118734329852813977L;

	public ReferenceBasePredictorNodeOnlyInterval29() {
	}

	@Override
	public CacheStrategySpec getCacheStrategy() {

		CacheStrategySpec css = new CacheStrategySpec(CacheStrategy.DENSE);		

		return css;
	}

}
