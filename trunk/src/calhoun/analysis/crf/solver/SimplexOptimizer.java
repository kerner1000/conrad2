package calhoun.analysis.crf.solver;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.CRFObjectiveFunctionGradient;
import calhoun.analysis.crf.CRFTraining;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.util.ErrorException;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import flanagan.math.Minimisation;
import flanagan.math.MinimisationFunction;

/** uses a nelder-mead algorithm (the simplex method) to do a general function optimization objective function.
 * This optimization method does not use the gradient but requires many iterations and so is mainly useful in
 * debugging.  Use the {@link StandardOptimizer} for most problems.
 * <p>
 * This optimizer has several configuration properties that allow control over the optimization process:
 * <ul>
 * <li> <b><code>maxIters</code></b> - The maximum number of iterations (objective function evaluations) to attempt
 * <li> <b><code>requireConverge</code></b> - if true, throws an error if convergence is not reached in the maximum number 
 * of iterations.  Otherwise, the current feature weights are returned when <code>maxIters</code> is reached.
 * <li> <b><code>stepSize</code></b> - the size of the initial changes made to the weights when exploring the objective function.
 * <li> <b><code>starts</code></b> - an initial set of guesses at feature weights.  Defaults to 1.0
 * </ul>
 * */
public class SimplexOptimizer implements CRFTraining {
	@SuppressWarnings("unused")
	private static final Log log = LogFactory.getLog(SimplexOptimizer.class);
	
	// Configuration
	CRFObjectiveFunctionGradient gradFunc;
	int maxIters = 500;
	double stepSize = 0.5;
	boolean requireConvergence = true;
	double[] starts = null;
	
	public double[] optimize(ModelManager fm, List<? extends TrainingSequence<?>> data) {
		gradFunc.setTrainingData(fm, data);		
		final int nFeatures = fm.getNumFeatures();
		MinimisationFunction mFunc = new MinimisationFunction() {
			double[] grad = new double[nFeatures];
			
			public double function(double[] d) {
				return -gradFunc.apply(d, grad);
			}
		};
		
		Minimisation m = new Minimisation();
		DenseDoubleMatrix1D steps = new DenseDoubleMatrix1D(nFeatures);
		m.setNmax(maxIters);
		if(starts == null) {
			starts = new double[nFeatures];
			Arrays.fill(starts, 1.0);
		}
		steps.assign(stepSize);
		m.nelderMead(mFunc, starts, steps.toArray());
		if(requireConvergence && !m.getConvStatus()) {
			throw new ErrorException("Convergence not reached.");
		}
		//m.print("test/working/nelder.txt");
		return m.getParamValues();
	}

	/** returns the configured objective function gradient which will be 
	 * used by the optimizer during the training process.
	 * @return the configured objective function gradient
	 */
	public CRFObjectiveFunctionGradient getObjectiveFunction() {
		return gradFunc;
	}
	
	/** sets the objective function gradient.  Called automatically during configuration. */
	public void setObjectiveFunction(CRFObjectiveFunctionGradient objectiveFunction) {
		this.gradFunc = objectiveFunction;
	}

	public int getMaxIters() {
		return maxIters;
	}

	public void setMaxIters(int maxIters) {
		this.maxIters = maxIters;
	}

	public boolean isRequireConvergence() {
		return requireConvergence;
	}

	public void setRequireConvergence(boolean requireConvergence) {
		this.requireConvergence = requireConvergence;
	}

	public double[] getStarts() {
		return starts;
	}

	public void setStarts(double[] starts) {
		this.starts = starts;
	}

	public double getStepSize() {
		return stepSize;
	}

	public void setStepSize(double stepSize) {
		this.stepSize = stepSize;
	}
}
