package calhoun.analysis.crf.test;

import java.io.IOException;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.Conrad;
import calhoun.util.AbstractTestCase;
import calhoun.util.Assert;

public class Interval13BaselineTest extends AbstractTestCase {
	private static final Log log = LogFactory.getLog(Interval13BaselineTest.class);

	/* ****************************************************
	 * These first tests show that various configurations of the interval13 model, both the real
	 * semi-Markov models and their toy versions that are Markov models, will run to completion through both training and testing.
	 * Making no effort yet to verify the correctness of results.
	 * 
	 * These first few tests are not very strong.
	 */

	public void testUntiedModel() throws Exception {
		// Tests Conrad on a SEMI-MARKOV model invoking it through COMMAND LINE
		Conrad.main(new String[] {"train", "test/input/interval13/config/baseline_untied.xml", "test/input/interval13/data/oneGeneTrain.interval13.txt", "test/working/interval13BaselineUntied.ser"});

		Conrad c = Conrad.read("test/working/interval13BaselineUntied.ser");
		double[] weights = c.getWeights();
		for(int i=0; i<weights.length; ++i) {
			log.info(c.getFeatureName(i)+": "+weights[i]);
		}
		assertEquals(17, weights.length);
	}
		
	public void testSemiMarkovCommandLine() throws Exception {
		// Tests Conrad on a SEMI-MARKOV model invoking it through COMMAND LINE

		Conrad.main(new String[] {"train", "test/input/interval13/config/semiMarkovZeroPad.xml", "test/input/interval13/data/oneGeneTrain.interval13.txt", "test/working/interval13BaselineModelTest2.ser"});

		Conrad.main(new String[] {"test", "test/working/interval13BaselineModelTest2.ser", "test/input/interval13/data/shortTrain.interval13.txt", "test/working/interval13BaselineModelTestPredicted.txt"});
	}
		
	public void testInference() throws Exception {
		
		// Do first using Semi-Markov model
		Conrad.main(new String[] {"train", "test/input/interval13/config/semiMarkovZeroPadNoTrain.xml", "test/input/interval13/data/oneGeneTrain.interval13.txt", "test/working/interval13BaselineModelTest.ser"});
		Conrad.main(new String[] {"test", "test/working/interval13BaselineModelTest.ser", "test/input/interval13/data/oneGeneTrain.interval13.txt", "test/working/interval13SemiMarkovModelTestPredicted.txt"} );
		
		
		// Do next using Markov model
		Conrad.main(new String[] {"train", "test/input/interval13/config/markovNoTrain.xml", "test/input/interval13/data/oneGeneTrain.interval13.txt", "test/working/interval13MarkovModelTest.ser"});
		Conrad.main(new String[] {"test", "test/working/interval13MarkovModelTest.ser", "test/input/interval13/data/oneGeneTrain.interval13.txt", "test/working/interval13MarkovModelTestPredicted.txt"} );
		
		assertFilesMatch("test/working/interval13SemiMarkovModelTestPredicted.txt","test/working/interval13MarkovModelTestPredicted.txt");	
	}

	public void testTrainingAndInferenceMarkovVsSemiMarkov() throws Exception {	
		String fileModel1 = "test/input/interval13/config/markov.xml";
		String fileModel2 = "test/input/interval13/config/semiMarkovZeroPad.xml";
		String fileData = "test/input/interval13/data/oneGeneTrain.interval13.txt";
			
		compareTwoEquivalentModelsTrainTestSameData(fileModel1,fileModel2,fileData);
	}

	public void testTrainingAndInferenceCacheProcessorCompare() throws Exception {
		/* SemiMarkovModel (with no length dependent features, just one NodeBoundary) 
		 * using Semi-Markov training and inference, using two different
		 * Cache Processors: AllSparse and CacheProcessorDeluxe
		 */ 
		String fileModel1 = "test/input/interval13/config/semiMarkovZeroPad.xml";
		String fileModel2 = "test/input/interval13/config/semiMarkovZeroPadCPD.xml";
		String fileData = "test/input/interval13/data/oneGeneTrain.interval13.txt";
			
		compareTwoEquivalentModelsTrainTestSameData(fileModel1,fileModel2,fileData);		
	}

	public void testMarkovModelWithMarkovTrainingTestingVsSemiMarkovCPDTrainingTesting() throws Exception {
		/* A markov model using simple engines and using semi-Markov training and CacheProcessorDeluxe
		 */
		String fileModel1 = "test/input/interval13/config/strictMarkovCPD.xml";
		String fileModel2 = "test/input/interval13/config/strictMarkov.xml";
		String fileData = "test/input/interval13/data/oneGeneTrain.interval13.txt";
			
		compareTwoEquivalentModelsTrainTestSameData(fileModel1,fileModel2,fileData);		
	}
	
	public void test2() throws Exception {
		/* Nontrivial test comparing Markov and Semi-Markov training,
		 * using same model and CacheProcessor and data.
		 */
		String fileModel1 = "test/input/interval13/config/strictMarkovCPD.xml";
		String fileModel2 = "test/input/interval13/config/strictMarkovCPD2.xml";
		String fileData = "test/input/interval13/data/oneGeneTrain.interval13.txt";
			
		compareTwoEquivalentModelsTrainTestSameData(fileModel1,fileModel2,fileData);		
	}
	
	public void test3() throws Exception {
		/* Compare two cache processors , CPD and Allsparse, in the case that
		 * the model includes a feature that is genuinely length dependent.
		 */
		String fileModel1 = "test/input/interval13/config/lengthDependentCPD.xml";
		String fileModel2 = "test/input/interval13/config/lengthDependentAllSparseNoMin.xml"; // Note: this model is very slow when building cache
		String fileData = "test/input/interval13/data/oneGeneTrain.interval13.txt";
			
		// Change the tolerance for the weights.  They can be different due to numerical error.
		compareTwoEquivalentModelsTrainTestSameData(fileModel1,fileModel2,fileData, 0.1);		
		
		
		// The block below will failif run because AllSparseLengthCacheProcessor does not yet support specifying the minimum statelength.
		/* Dave needs to fix this bc probably a mistake in AllSparseCacheProcessor not looking at min-lengths.
		 * Only difference from test3 is in the ocnfig file for CPD I uncomment the minlengths.
		 */
		//String fileModel1a = "test/input/interval13/config/lengthDependentCPDuncommented.xml";
		//String fileModel2a = "test/input/interval13/config/lengthDependentAllSparse.xml"; // Note: this model is very slow when building cache
		//String fileDataa = "test/input/interval13/data/oneGeneTrain.interval13.txt";
			
		//compareTwoEquivalentModelsTrainTestSameData(fileModel1a,fileModel2a,fileDataa, 0.1);		
		
		
	}
	
	public void testReadingBaselineXMLConfigFile() throws Exception {
		Conrad cr = new Conrad("test/input/interval13/config/ssbaselineSplitInput.xml");
		System.out.println("Trying to train model");
		cr.train("test/input/interval13/data/splitInputOneGeneTrain"); 
		System.out.println("Trying to test model");
		cr.test("test/input/inputFilesTest/");
	}

	public void testReadingSSBaselineCPDXMLConfigFile() throws Exception {
		Conrad cr = new Conrad("test/input/interval13/config/ssbaselineCPD.xml");
		System.out.println("Trying to train model");
		cr.train("test/input/interval13/data/oneGeneTrain.interval13.txt");
		System.out.println("Trying to test model");
		cr.test("test/input/interval13/data/oneGeneTrain.interval13.txt");
	}
	
	public void compareTwoEquivalentModelsTrainTestSameData( String fileModel1, String fileModel2, String fileData) throws IOException {
		compareTwoEquivalentModelsTrainTestSameData(fileModel1, fileModel2, fileData, 0.001);
	}
	
	public void compareTwoEquivalentModelsTrainTestSameData( String fileModel1, String fileModel2, String fileData, double tolerance ) throws IOException {
		Conrad cr1 = new Conrad(fileModel1);
		cr1.train(fileData);
		double[] weights1 = cr1.getWeights();	
		
		Conrad cr2 = new Conrad(fileModel2);
		cr2.train(fileData);
		double[] weights2 = cr2.getWeights();	
		
		Assert.a(weights2.length == weights1.length,"length1 is " + weights1.length  + "   and weights2 is " + weights2.length);
		assertArrayEquals(weights2,weights1, tolerance);		

		cr1.test(fileData);	
		cr2.test(fileData);				

		String s1 = cr1.getOutputHandler().toString();
		System.out.println("String of output from result1 is");
		System.out.println(s1);
		
		String s2 = cr2.getOutputHandler().toString();
		System.out.println("String of output from result2 is");
		System.out.println(s2);
		
		assertEquals(s1,s2);
	}
	
}
