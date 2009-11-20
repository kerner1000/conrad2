package calhoun.analysis.crf;

import calhoun.analysis.crf.io.InputSequence;

/** holds features whose value depends on the current state and the previous state.
 * This is the most general type of feature for a markov CRF.
 */
public interface FeatureManagerEdge<InputType> extends FeatureManager<InputType> {

	/** Evaluates the set of features managed by this object for the given arguments. */
	void evaluateEdge(InputSequence<? extends InputType> seq, int pos, int prevState, int state, FeatureList result);
}
