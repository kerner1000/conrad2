package calhoun.analysis.crf;

import calhoun.analysis.crf.io.InputSequence;


/** holds features whose value only depends on the current state, and not the previous state.
 * An {@link FeatureManagerEdge} is more general than a node feature, but node features are more efficient
 * since the solver does not need to evaluate them separately for each transition.  
 */
public interface FeatureManagerNode<InputType> extends FeatureManager<InputType> {

	/** Evaluates the set of features managed by this object for the given arguments. */
	void evaluateNode(InputSequence<? extends InputType> seq, int pos, int state, FeatureList result);
}
