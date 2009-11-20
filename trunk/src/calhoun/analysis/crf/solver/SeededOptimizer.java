package calhoun.analysis.crf.solver;

import java.io.File;
import java.util.List;

import calhoun.analysis.crf.CRFTraining;
import calhoun.analysis.crf.Conrad;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.io.TrainingSequence;

/** an optimizer that uses the weights from another model as the seed for a new optimization.  Allows a second pass optimization on an already trained model. */ 
public class SeededOptimizer implements CRFTraining {

	CRFTraining seededOptimizer;
	File seedModel;

	public double[] optimize(ModelManager fm, List<? extends TrainingSequence<?>> data) {
		try {
			Conrad seed = Conrad.read(seedModel.getPath());
			seededOptimizer.setStarts(seed.getWeights());
			return seededOptimizer.optimize(fm ,data);
		}
		catch(Exception ex) {
			throw new RuntimeException(ex);
		}
	}

	/** starting weights are ignored for a SeededOptimizer, since they are always taken from the seeded model. */  
	public void setStarts(double[] weights) {
		throw new UnsupportedOperationException("Cannot set starts on a SeededOptimizer.");
	}
	
	/**
	 * @return the seededOptimizer
	 */
	public CRFTraining getSeededOptimizer() {
		return seededOptimizer;
	}

	/**
	 * @param seededOptimizer the seededOptimizer to set
	 */
	public void setSeededOptimizer(CRFTraining seededOptimizer) {
		this.seededOptimizer = seededOptimizer;
	}

	/**
	 * @return the seedModel
	 */
	public File getSeedModel() {
		return seedModel;
	}

	/**
	 * @param seedModel the seedModel to set
	 */
	public void setSeedModel(File seedModel) {
		this.seedModel = seedModel;
	}
}
