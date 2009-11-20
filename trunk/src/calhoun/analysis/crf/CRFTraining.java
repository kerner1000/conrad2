package calhoun.analysis.crf;

import java.util.List;

import calhoun.analysis.crf.io.TrainingSequence;

/** an interface to numerical solvers for optimizing the CRF objective function.  Conrad will call this when training the weights,
 * and this class is reponsible for evlauating the objective function iteratively and determining the weights.  Usually,
 * this class will just be a wrapper around some standard numerical solving package.
 * <p>
 * The optimizer is not required to find an optimal set of weights.  Although this is usually the goal and is usually feasible, the
 * interface only requires that a feature weight be assigned for each feature.  It is legal for the optimizer to return suboptimal weights.
 */
public interface CRFTraining {

	/** find the set of weights which maximizes the value of the objective function.  The {@link CRFObjectiveFunctionGradient} object will already
	 * be configured witht eh appropriate training data, so it appears here as a pure function evaluation.
	 *   
	 * @param fm		the model to train on.
	 * @param data		the training data to use for training.
	 * @return an array of doubles containing the optimal weights.
	 */
	double[] optimize(ModelManager fm, List<? extends TrainingSequence<?>> data);

	/** Sets the starting weights for the optimization. 
	 * @param weights an array of weights to use as the starting point for the optimizer.  The optimzier is not required to use this as a starting point.*/ 
	void setStarts(double[] weights);
}