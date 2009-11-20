package calhoun.analysis.crf.solver;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import riso.numerical.LBFGS;
import calhoun.analysis.crf.CRFObjectiveFunctionGradient;
import calhoun.analysis.crf.CRFTraining;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.util.Assert;
import calhoun.util.ErrorException;

/** uses a L-BFGS algorithm to optimize the objective function.  This is the algorithm normally used. 
 * Several configuration parameters all you to control optimization parameters:
 * 
 * <p>
 * This optimizer requires one property to be set:
 * <p>
 * <b><code>gradFunc</code></b> - This is the gradient function that the optimizer will use.  It must be a reference to a bean that implements the {@link CRFObjectiveFunctionGradient} interface. 
 * <p>
 * In addition, it has several optional properties that allow control over the optimization process:
 * <ul>
 * <li> <b><code>debugLevel</code></b> - explicitly sets the debug level in the underlying L-BFGS solver.  Normally it is set based
 * on the debug level for this object. 
 * <li> <b><code>epsForConvergence</code></b> - epsilon value used to determine when to halt the optimzation and declase convergence.  Defaults
 * to 0.0001.
 * <li> <b><code>fixFirstWeight</code></b> - if true, the first feature weight will be fixed at 1.0 and will not be allowed to change 
 * during the optimization
 * <li> <b><code>maxIters</code></b> - The maximum number of iterations (objective function evaluations) to attempt
 * <li> <b><code>mForHessian</code></b> - sets the <code>mForHession</code> parameter in the underlying L-BFGS solver.  Defaults to 20.
 * <li> <b><code>quadraticRegularization</code></b> - if set to a nonzero value, regularizes feature weights by imposing 
 * a penalty on the objective function based on the absolute sizes of the weights.  Can help in cases where weights
 * get very large.
 * <li> <b><code>requireConverge</code></b> - if true, throws an error if convergence is not reached in the maximum number 
 * of iterations.  Otherwise, the current feature weights are returned when <code>maxIters</code> is reached.
 * <li> <b><code>starts</code></b> - an initial set of guesses at feature weights.  Defaults to 1.0
 * <li> <b><code>unchangedObjective</code></b> - the number of times the objectiveFunction must return an identical value before the optimization is assumed to have converged.  
 *   This protects against a situation where the epsForConvergence is set low enough that even the smallest possible changes in weights can't find an optimum.  Defaults to 5
 * </ul>
 * */
public class StandardOptimizer implements CRFTraining {
	private static final Log log = LogFactory.getLog(StandardOptimizer.class);
	boolean debug = log.isDebugEnabled();
	
	// Configuration
	CRFObjectiveFunctionGradient gradFunc;
	int maxIters = 2000;
	int mForHessian = 20;
	int debugLevel = 0;
	boolean requireConvergence = true;
	double epsForConvergence = 0.00001;
	double[] starts = null;
	boolean fixFirstWeight = false;
	double quadraticRegularization = 0.0;
	int unchangedObjective = 5;
	
	public double[] optimize(ModelManager fm, List<? extends TrainingSequence<?>> data) {
		gradFunc.setTrainingData(fm, data);
		int nFeatures = fm.getNumFeatures();
		if (requireConvergence) {
			log.info("NOTE: You ARE requiring convergence of LBFGS");
		} else {
			log.warn("NOTE: You ARE NOT requiring convergence of LBFGS");			
		}

		double f = Float.NaN, xtol = 1.0e-16; // machine precision
		int iprint[] = new int[2];
		iprint[0] = debugLevel - 2;
		iprint[1] = debugLevel - 1;
		int[] iflag = new int[] {0};
		int numFeatures = nFeatures;

		double[] lambda = new double[nFeatures];
		double[] grad = new double[nFeatures];

		if(fixFirstWeight) {
			log.info("Fixing first weight to 1.0.");
			numFeatures -= 1;
		}
		
		double[] optLambda = new double[numFeatures];
		double[] optGrad = new double[numFeatures];
		double[] diag = new double[numFeatures]; // needed by the optimizer

		if(starts == null) {
			Arrays.fill(lambda, 1.0);
		}
		else {
			Assert.a(starts.length == lambda.length, "Received ", starts.length, " initial weights.  Expected: ", lambda.length);
			System.arraycopy(starts, 0, lambda, 0, lambda.length);
		}

		int icall = 0;
		float lastObjective = Float.NaN;
		int runningObjective = 0;
		do {
			if(fixFirstWeight) {
				lambda[0] = 1.0;
			}
			
			try {
				f = gradFunc.apply(lambda, grad);
			}
			catch(RuntimeException ex) {
				if(requireConvergence) {
					throw ex;
				} else {
					log.warn("Exception thrown while calculating gradient.  Possible numeric problem.");
					log.warn(ex);
					return lambda;
				}
			}

			// Take the gradient, normalize by the total length of the sequence, and make it a minimization instead of a maximization problem.
			//f = -f; 
			f = -f; 
			if(fixFirstWeight) {
				for (int j = 1; j < lambda.length; j++) {
					optGrad[j-1] = -grad[j];
					//optGrad[j-1] = -grad[j];
				}
				for (int j = 1; j < lambda.length; j++) {
					optLambda[j-1] = lambda[j];
				}
			}
			else {
				for (int j = 0; j < lambda.length; j++) {
					optGrad[j] = -grad[j];
					//optGrad[j] = -grad[j];
				}
				for (int j = 0; j < lambda.length; j++) {
					optLambda[j] = lambda[j];
				}
			}
			// Add a regularization term quadraticRegularization*sum(lambda_i^2)
			for (int j=0; j<lambda.length; j++) {
				f += quadraticRegularization*lambda[j]*lambda[j];
				optGrad[j] += 2*quadraticRegularization*lambda[j];
			}
			
			// Check for an objective value that hasn't changed over several iterations.
			if(((float) f) == lastObjective) {
				++runningObjective;
				log.info("Objective value unchanged: "+lastObjective+" returned "+runningObjective+" times.");
				if(runningObjective > 0 && runningObjective >= unchangedObjective) {
					log.warn("Same objective value: "+lastObjective+" returned "+(unchangedObjective+1)+" times.  Assuming convergence.");
					return lambda;
				}
			}
			else {
				runningObjective = 0;
				lastObjective = (float) f;
			}
			
			try {
				LBFGS.gtol = 0.1;
				LBFGS.lbfgs(numFeatures, mForHessian, optLambda, f, optGrad, false, diag, iprint, epsForConvergence, xtol, iflag);
			} catch (LBFGS.ExceptionWithIflag e) {
				if(requireConvergence)
					throw new ErrorException("lbfgs failed", e);
				else {
					log.warn("Lbfgs failed, proceeding anyway", e);
					break;
				}
			}
			if(fixFirstWeight) {
				for (int j = 1; j < lambda.length; j++) {
					lambda[j] = optLambda[j-1];
				}
			}
			else {
				for (int j = 0; j < lambda.length; j++) {
					lambda[j] = optLambda[j];
				}
			}
			icall += 1;
		} while ((iflag[0] != 0) && (icall < maxIters));
		if(requireConvergence && !(iflag[0] == 0)) {
			throw new ErrorException("Convergence not reached.");
		}
		return lambda;
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

	public int getDebugLevel() {
		return debugLevel;
	}

	public void setDebugLevel(int debugLevel) {
		this.debugLevel = debugLevel;
	}

	public double getEpsForConvergence() {
		return epsForConvergence;
	}

	public void setEpsForConvergence(double epsForConvergence) {
		this.epsForConvergence = epsForConvergence;
	}

	public boolean isFixFirstWeight() {
		return fixFirstWeight;
	}

	public void setFixFirstWeight(boolean fixFirstWeight) {
		this.fixFirstWeight = fixFirstWeight;
	}

	public int getMaxIters() {
		return maxIters;
	}

	public void setMaxIters(int maxIters) {
		this.maxIters = maxIters;
	}

	public int getMForHessian() {
		return mForHessian;
	}

	public void setMForHessian(int forHessian) {
		mForHessian = forHessian;
	}

	public double getQuadraticRegularization() {
		return quadraticRegularization;
	}

	public void setQuadraticRegularization(double quadraticRegularization) {
		this.quadraticRegularization = quadraticRegularization;
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

	/**
	 * @return Returns the unchangedObjective.
	 */
	public int getUnchangedObjective() {
		return unchangedObjective;
	}

	/**
	 * @param unchangedObjective The unchangedObjective to set.
	 */
	public void setUnchangedObjective(int unchangedObjective) {
		this.unchangedObjective = unchangedObjective;
	}
}
