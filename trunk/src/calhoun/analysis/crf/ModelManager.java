package calhoun.analysis.crf;

import calhoun.util.DenseBooleanMatrix2D;

/** overall itnerface for an entire model.  The model contains the hidden state defintions and transitions as well as all of the features.
 */ 
public interface ModelManager extends FeatureManagerEdge, FeatureManagerNode, FeatureManagerEdgeExplicitLength, FeatureManagerNodeExplicitLength {
	
	/** the number of hidden states.
	 * @return the number of hidden states
	 */
	public int getNumStates();

	/** looks up the human readable name of a hidden state given it's index.
	 * @return string name of the state with this index.
	 */
	public String getStateName(int state);

	/** looks up a hidden state's index given it's human readable name.
	 * @return index of the state with this name.
	 */
	public int getStateIndex(String name);
	
	/** returns a boolean matrix where each <code>(row, column)</code> entry contains true if the
	 * model has a legal transition between state <code>row</code> and state <code>column</code>.
	 * A ModelManager may safely return null if all transitions are legal.
	 * @return a boolean matrix containing the legal transitions.
	 */
	public DenseBooleanMatrix2D getLegalTransitions();
}
