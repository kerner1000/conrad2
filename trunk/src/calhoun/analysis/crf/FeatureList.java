package calhoun.analysis.crf;

/** used to pass feature evaluations from a FeatureManager to the Conrad engine during an <code>evaluate</code> call.
 * A <code>FeatureList</code> is always the last argument passed to an <code>evaluate</code> call on <code>FeatureManager</code>.
 * It allows the feature to add in evaluations or declare the given node or edge invalid.  This interface is present so that
 * new memory does not have to be allocated during the evaluation process and so that a given FeatureManager can return any
 * arbitrary number of feature evaluations from an <code>evaluate</code> call.
 * 
*/
public interface FeatureList {
	/** returns a feature value for this evaluation.  Since a <code>FeatureManager</code> can provide evaluations for more than one
	 * feature, it is necessary to identify the specific feature by passing its <code>index</code>.<p>
	 * If an <code>evaluate</code> call is made and <code>addFeature</code> is not called for a given feature index, then the value for
	 * that feature is assumed to be zero.  This makes it very efficient to evaluate sparse features.<p> 
	 * The value can be any real number, and passing zero is allowed but unnecessary.
	 * @param index the index of the feature which has been evaluated
	 * @param value the value of the feature.
	 */ 
	void addFeature(int index, double value);
	
	/** specifies that node, edge, or segment specified in the <code>evaluate</code> call is invalid.  When invalid nodes or edges are found
	 * the engine will exclude all hidden state paths including that node, edge, or segment from both training and inference.  This is very
	 * different from returning zero for a given feature.  Note that if even one <code>FeatureManager</code> invalidates a node or edge, then
	 * no feature evaluations at that position will be used.  This is very different from returning zero, which simply means that the path is 
	 * valid, but the feature doesn't provide any information.
	 */
	void invalidate();

	/** checks if the current node, edge, or segment is valid.  This can be used to speed up processing in more complicated FeatureManagers, 
	 * since once invalidate has been called during an <code>evaluate</code> function there is no point in adding any new feature values.  
	 * @return true if the <code>FeatureList</code> is still valid.
	 */
	boolean isValid();
}
