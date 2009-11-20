package calhoun.analysis.crf;


/** An extension of FeatureManagerNode for situations where the nodes use in semi-markov models is non-standard.
 * It allows an arbitrary region at the right and left edge of each segment to be eliminated, presumable because a
 * different probabilistic model will be used at that position.
 */
public interface FeatureManagerNodeBoundaries<InputType> extends FeatureManagerNode<InputType> {
	/** A NodeBoundary defines how a node feature should be used in a semi-markov model.  The default behavior is that
	 * the value of a node will be the sum of the values of the node feature at all positions in the segment.  Using
	 * nodeBoundaries allows more flexibility.  
	 * This is somewhat ill-defined right now.  Designed specifically for EmissonMarkovFeatures in the semi-markov gene caller.*/ 
	
	/*
	 *  NOTE: These details below have now been wrapped up in the details of CacheStrategySpec, but we still want to keep this as a separate class.
	public static class NodeBoundary {
		int featureIndex;
		int state;
		int lookupTable;
		int leftBoundary;
		int rightBoundary;
		
		public NodeBoundary(int featureIndex, int state, int lookupTable, int leftBoundary, int rightBoundary) {
			featureIndex = this.featureIndex;
			state = this.state;
			lookupTable = this.lookupTable;
			leftBoundary = this.leftBoundary;
			rightBoundary = this.rightBoundary;
		}
	}
	

	List<NodeBoundary> getBoundaries();
	
	*/
}
