package calhoun.analysis.crf.solver;

import java.util.List;

import calhoun.analysis.crf.CRFTraining;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.io.TrainingSequence;

/** a two pass optimizer that does an initial optimization and then uses the weights generated from that as the start of a second optimization.  It is configured with
 * two separate CRFOptimizer oibjects for the two passes and additionally needs a CRFObjectiveFunctionGradient for the first pass. 
 * <p>
 * This optimizer has two configuration properties that must be set and control the optimization process
 * <ul>
 * <li> <b><code>firstPass</code></b> - A reference to a configured {@link CRFTraining} object that will be used to optimize the 
 * weights in the first pass of the optimization. 
 * <li> <b><code>secondPass</code></b> - A reference to a configured {@link CRFTraining} object that will be used to optimize the 
 * weights in the second pass of the optimization.  The second pass will begin using the weights computed during the first
 *  pass.
 */
public class TwoPassOptimizer implements CRFTraining {

	CRFTraining firstPass;

	CRFTraining secondPass;

	public double[] optimize(ModelManager fm, List<? extends TrainingSequence<?>> data) {
		double[] startingWeights = firstPass.optimize(fm, data);

		secondPass.setStarts(startingWeights);
		return secondPass.optimize(fm, data);
	}

	public void setStarts(double[] weights) {
		firstPass.setStarts(weights);
	}

	/**
	 * @return the firstPass
	 */
	public CRFTraining getFirstPass() {
		return firstPass;
	}

	/**
	 * @param firstPass the firstPass to set
	 */
	public void setFirstPass(CRFTraining firstPass) {
		this.firstPass = firstPass;
	}

	/**
	 * @return the secondPass
	 */
	public CRFTraining getSecondPass() {
		return secondPass;
	}

	/**
	 * @param secondPass the secondPass to set
	 */
	public void setSecondPass(CRFTraining secondPass) {
		this.secondPass = secondPass;
	}
}
