package calhoun.analysis.crf;

import calhoun.analysis.crf.CacheStrategySpec.CacheStrategy;

/** base class for feature implementations.  Takes care of some of the details and bookkeeping involved in
 * writing a feature.
 */
public abstract class AbstractFeatureManager<InputType> implements FeatureManager<InputType> {

	String inputComponent;

	public CacheStrategySpec getCacheStrategy() {
		return new CacheStrategySpec(CacheStrategy.UNSPECIFIED);
	}

	public String getInputComponent() {
		return inputComponent;
	}

	public void setInputComponent(String inputComponent) {
		this.inputComponent = inputComponent;
	}
}
