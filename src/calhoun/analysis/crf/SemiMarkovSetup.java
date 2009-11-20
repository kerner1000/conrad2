package calhoun.analysis.crf;

/** holds additional configuration information used in semi-Markov CRFs.  
 */
public class SemiMarkovSetup {
	boolean ignoreSemiMarkovSelfTransitions;
	private short[] maxLengths;
	private short[] minLengths;

	/** default constructor used by Spring during automatic configuration */
	public SemiMarkovSetup() { }
	/** constructs using this set of maximum lengths. */
	public SemiMarkovSetup(short[] lengths) { this.maxLengths = lengths; }
	/** constructs using this set of maximum lengths and value for the ignore flag. */
	public SemiMarkovSetup(short[] lengths, boolean ignore) { 
		this(lengths);
		ignoreSemiMarkovSelfTransitions = ignore;
	}
	
	/** constructs will all parameters. */
	public SemiMarkovSetup(short[] minLen, short maxLen[], boolean ignore) {
		this.minLengths = minLen;
		this.maxLengths = maxLen;
		this.ignoreSemiMarkovSelfTransitions = ignore;
	}
	
	/** gets the minimum lengths for each state.  The minimum value is 1.  This array defaults to
	 * all 1's if it was not set explicitly.
	 * @return an array of the minimum lengths of each state.
	 */  
	public short[] getMinLengths() {
		// Default minLengths to all zeros (and the same size as maxLengths).
		if(minLengths == null)
			minLengths = new short[maxLengths.length];
		return minLengths;
	}

	/** sets the minimum lengths for each state.
	 * @param lengths an array of the minimum lengths of each state.
	 */  
	public void setMinLengths(short[] lengths) {
		this.minLengths = lengths;
	}

	/** gets the maximum lengths for each state.  The minimum value is 1, which corresponds to the 
	 * special case of a Markov feature.  This value has no default.  It is usually configured through the XML model file.
	 * @return an array of the maximum lengths of each state.
	 */  
	public short[] getMaxLengths() {
		return maxLengths;
	}

	/** sets the maximum lengths for each state.
	 * @param lengths an array of the minimum lengths of each state.
	 */  
	public void setMaxLengths(short[] lengths) {
		this.maxLengths = lengths;
	}

	/** gets the flag indicating that self-transition edges for the semi-Markov states should be ignored.  This flag is useful
	 * if you are using the same model definition for both Markov and semi-Markov models.  In a Markov model a self-transition 
	 * means that a state can be repeated multiple times, essentially forming a segment.  However, in a semi-Markov a self-transition 
	 * means a transition between two segments with the same state.  This distinction means that the same state transition graph
	 * will mean different things under the two models.  If this flag is set to true, which is the 
	 * default, then it is assumed that no two segments in the semi-Markov model can have the same state, and so self-transitions in
	 * the semi-Markov states are ignored.  This forces the state transition graph to have the same meaning for both cases.
	 * <p>
	 * This flag does not affect Markov states, those with a max length of 1.  This means that a state may legally occur in consecutive
	 * positions if the maximum length for that state is one and the graph contains a self-transition for that state.
	 * @return true if self-transition edges in the model should be ignored.
	 */  
	public boolean isIgnoreSemiMarkovSelfTransitions() {
		return ignoreSemiMarkovSelfTransitions;
	}

	/** sets the flag indicating that self-transition edges for the semi-Markov states should be ignored.
	 * @param ignoreSemiMarkovSelfTransitions true if self-transition edges for the semi-Markov states should be ignored.
	 */  
	public void setIgnoreSemiMarkovSelfTransitions(
			boolean ignoreSemiMarkovSelfTransitions) {
		this.ignoreSemiMarkovSelfTransitions = ignoreSemiMarkovSelfTransitions;
	}
}
