/**
 * 
 */
package calhoun.analysis.crf;

import calhoun.analysis.crf.io.InputSequence;

/** an interface to inference algorithms for CRFs.  Given a model, a set of feature weights, and a set of input data, the
 * algorithm selects a sequence of hidden states.  The Viterbi dynamic programming algorithm and its variants are usually used
 * for this problem. 
 */
public interface CRFInference {
	
	/** holder which contains the results of an inference run.  The indexes of the predict hidden states are stored in the 
	 * <code>hiddenStates</code> array.  The best scores array is a column major array of the best scores to each estate and position.
	 */
	public static class InferenceResult {
		public int[] hiddenStates;
		public double[] bestScores;
	}
	
	/** Return the labelling that maximizes the conditional probability P(y|x).
	 * @param mm 		model to use for training  
	 * @param data 		input sequence to label
	 * @param weights	array of feature weights.  Usually these will be derived from a training pass.
	 * @return inference result containing the hidden states which are predicted and the score outputs.
	 */
	InferenceResult predict(ModelManager mm, InputSequence<?> data, double[] weights);
}