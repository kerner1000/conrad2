package calhoun.analysis.crf.solver;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.CRFTraining;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.util.Assert;

/** a dummy optimizer that just fixed the weights at values specified in the configuration. 
 * If no weights are specified, then a default all weights being fixed at 1.0 is used. */
public class FixedWeightOptimizer implements CRFTraining {
	private static final Log log = LogFactory.getLog(FixedWeightOptimizer.class);

	double[] starts;
	
	public double[] optimize(ModelManager fm, List<? extends TrainingSequence<?>> data) {
		int N = fm.getNumFeatures();
		if (starts == null) {
			log.warn("FixedWeightOptimizer called but weights not specified; setting all weights to 1.0");
			starts = new double[N];
			Arrays.fill(starts, 1.0);
		}
		Assert.a(starts.length == N , "Number of specified starting weights must equal the number of Features");
		return starts;
	}

	/** gets the fixed weights that will be used in place of an optimization.
	 * @return array of feature weights to use */
	public double[] getStarts() {
		return starts;
	}

	/** sets the values to use as feature weights.  A null value sets all weights to 1.0.
	 * @param starts array of feature weights to use */
	public void setStarts(double[] starts) {
		this.starts = starts;
	}
}
