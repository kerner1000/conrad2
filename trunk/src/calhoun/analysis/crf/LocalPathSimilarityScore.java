package calhoun.analysis.crf;

import calhoun.analysis.crf.io.TrainingSequence;

/** interface for local scoring functions used by the local score gradient functions.
 * The alternate objective functions that use a local similarity score to calculate use this
 * interface to access the local scoring functions.
 */
public interface LocalPathSimilarityScore {
	
	/** compute a real-valued score between a given path and the true path at a given position. 
	 * 
	 * @param yprev		the hidden state at the previous position in this path
	 * @param y			the hidden state at the current position in this path
	 * @param seq		the training sequence containing the input and the true path
	 * @param pos		the position at which the score is to be calculated
	 * @return the similarity score computed at this position for this path
	 */
	public double evaluate(int yprev, int y, TrainingSequence<?> seq, int pos);
}
