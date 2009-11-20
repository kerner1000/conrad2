package calhoun.analysis.crf;

import calhoun.analysis.crf.io.InputSequence;


/** Since Features are computed as needed dynamically at runtime, a FeatureManager controls training and dynamic creation of the features. 
 * A FeatureManager manages features generated from a given InputType. 
 * 
 * This is a feature manager that implements an explicit length distribution.
 * */
public interface FeatureManagerNodeExplicitLength<InputType> extends FeatureManager<InputType> {

	void evaluateNodeLength(InputSequence<? extends InputType> seq, int pos, int length, int state, FeatureList result);
}
