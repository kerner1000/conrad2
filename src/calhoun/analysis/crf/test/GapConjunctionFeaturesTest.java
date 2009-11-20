package calhoun.analysis.crf.test;


import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.Conrad;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.util.AbstractTestCase;

public class GapConjunctionFeaturesTest extends AbstractTestCase {
	private static final Log log = LogFactory.getLog(GapConjunctionFeaturesTest.class);
	boolean debug = log.isDebugEnabled();
	
	public void testGapFeature() throws Exception {

		String configFile = "test/input/multipleAlignmentUnitTest/basic_ModelDSComp.xml";
		String inputFile = "test/input/multipleAlignmentUnitTest/tiny_5way.txt";
		
		Conrad crf = new Conrad(configFile);
		List<? extends TrainingSequence<?>> t = crf.getInputHandler().readTrainingData(inputFile);

		crf.train(t);
		double[] weights = crf.getWeights();

		System.out.println("The weights are: ");
		for (int i=0; i<weights.length; i++ ) { System.out.print("\t"+weights[i]);		}
		System.out.println();
	}
	
}
