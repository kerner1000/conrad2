/**
 * 
 */
package calhoun.analysis.crf;

import java.util.List;

import calhoun.analysis.crf.io.TrainingSequence;

/** an interface to algorithms that compute an objective function and its gradient for CRFs.  The objective function is used
 * to determine the optimal weights during the CRF training process.  Usually the objective function value and gradient
 * will be computed many times for different values of feature weights during a training.  A numerical optimizer which implements
 * the {@link CRFTraining} interface is responsible for choosing the weights and calling this class.
 * 
 *  Different implementations are avialable to handle the Markov and semi-Markov CRFs and to handle different speed/memory tradeoffs.
 */
public interface CRFObjectiveFunctionGradient {

	/** sets the training data that will be used for evaluation of the objective function.  This function will be called before <code>apply</code>
	 * is called to set up the training data.  Since it is expected that <code>apply</code> will be called many times, this funtion is 
	 * the place to do one time setup and caching.
	 * 
	 * @param fm	the model to use.  Defines the hidden states, transitions, and features.
	 * @param data	the training sequences on which to calculate the objective function.
	 */
	void setTrainingData(ModelManager fm, List<? extends TrainingSequence<?>> data);

	/** computes the objective function value and the gradient.  This function will be called for each iteration of the numerical solver during training.
	 * In each iteration it takes a new set of feature weights and computes a value for the objective function along with its gradient.
	 * 
	 * @param weights an array of feature weights to use.  The length will equal the number of features in the model, and the values will change 
	 * for each call to apply. 
	 * @param grad an array which must be filled with the gradient vector when the function returns.  For each feature index, the array should contain an
	 * entry with the partial derivative with respect to that feature.
	 * @return the value of the objective function.  This is what the numerical optimizer will attempt to maximize.
	 */
	double apply(double[] weights, double[] grad);

	/** Frees resources allocated by setTrainingData */
	void clean();
}