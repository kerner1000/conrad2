package calhoun.analysis.crf.test;

import java.io.IOException;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.Conrad;
import calhoun.analysis.crf.solver.MaximumLikelihoodGradient;
import calhoun.analysis.crf.solver.MaximumLikelihoodSemiMarkovGradient;
import calhoun.analysis.crf.solver.StandardOptimizer;
import calhoun.util.AbstractTestCase;
import calhoun.util.Assert;

public class Interval29BaselineTest extends AbstractTestCase {
	private static final Log log = LogFactory.getLog(Interval29BaselineTest.class);
	public void testDummy() throws Exception {
		
	}
	public void testMarkov13VsSemiMarkov29() throws Exception {	
		String fileModel1 = "test/input/interval13/config/markov.xml";
		String fileModel2 = "test/input/interval29/config/semiMarkovZeroPad.xml";
		String fileData = "test/input/interval13/data/oneGeneTrain.interval13.txt";
			
		compareTwoEquivalentModelsTrainTestSameDataWithDifferentStates(fileModel1,fileModel2,fileData,0.01);
	}
	
	public void testMarkov13VsMarkov29() throws Exception {
		String fileModel1 = "test/input/interval13/config/strictMarkovCPD.xml";
		String fileModel2 = "test/input/interval29/config/strictMarkovCPDInt29.xml";
		String fileData = "test/input/interval13/data/oneGeneTrain.interval13.txt";
			
		compareTwoEquivalentModelsTrainTestSameDataJustFeatureSums(fileModel1,fileModel2,fileData,0.0001);		
	}

//public void testSemiMarkov13VsSemiMarkov29() throws Exception {
//	// This is a great test, but it's really slow and shouldn't give us any more 
//	// information than the Markov test above
//	String fileModel1 = "test/input/interval13/config/strictMarkovCPD2.xml";
//	String fileModel2 = "test/input/interval29/config/strictMarkovCPD2Int29.xml";
//	String fileData = "test/input/interval13/data/oneGeneTrain.interval13.txt";
//		
//	compareTwoEquivalentModelsTrainTestSameDataJustFeatureSums(fileModel1,fileModel2,fileData,0.001);		
//}
			
	public void compareTwoEquivalentModelsTrainTestSameDataJustFeatureSums( String fileModel1, String fileModel2, String fileData, double tolerance ) throws IOException {
		Conrad cr1 = new Conrad(fileModel1);
		cr1.train(fileData);
		double[] weights1 = cr1.getWeights();	
		
		double[] sums1;
		StandardOptimizer stdOpt = (StandardOptimizer)cr1.getOptimizer();
		if (stdOpt.getObjectiveFunction() instanceof MaximumLikelihoodGradient) {
			MaximumLikelihoodGradient grad = (MaximumLikelihoodGradient)stdOpt.getObjectiveFunction();
			sums1 = grad.getFeatureSums();
		} else {
			Assert.a(stdOpt.getObjectiveFunction() instanceof MaximumLikelihoodSemiMarkovGradient);
			MaximumLikelihoodSemiMarkovGradient grad = (MaximumLikelihoodSemiMarkovGradient)stdOpt.getObjectiveFunction();
			sums1 = grad.getFeatureSums();
		}
		
		cr1.test(fileData);	
		
		String s1 = cr1.getOutputHandler().toString();
		System.out.println("String of output from result1 is");
		System.out.println(s1);
		
		cr1 = null;
		
		Conrad cr2 = new Conrad(fileModel2);
		cr2.train(fileData);
		double[] weights2 = cr2.getWeights();	
		
		double[] sums2;
		stdOpt = (StandardOptimizer)cr2.getOptimizer();
		if (stdOpt.getObjectiveFunction() instanceof MaximumLikelihoodGradient) {
			MaximumLikelihoodGradient grad = (MaximumLikelihoodGradient)stdOpt.getObjectiveFunction();
			sums2 = grad.getFeatureSums();
		} else {
			Assert.a(stdOpt.getObjectiveFunction() instanceof MaximumLikelihoodSemiMarkovGradient);
			MaximumLikelihoodSemiMarkovGradient grad = (MaximumLikelihoodSemiMarkovGradient)stdOpt.getObjectiveFunction();
			sums2 = grad.getFeatureSums();
		}
		
		cr2.test(fileData);	
		
		String s2 = cr2.getOutputHandler().toString();
		System.out.println("String of output from result2 is");
		System.out.println(s2);
		
		Assert.a(weights2.length == weights1.length,"length1 is " + weights1.length  + "   and weights2 is " + weights2.length);
		
		assertArrayEquals(sums2,sums1, tolerance);
	}
	
	public void compareTwoEquivalentModelsTrainTestSameData( String fileModel1, String fileModel2, String fileData, double tolerance ) throws IOException {
		Conrad cr1 = new Conrad(fileModel1);
		cr1.train(fileData);
		double[] weights1 = cr1.getWeights();	
		
		cr1.test(fileData);	
		
		String s1 = cr1.getOutputHandler().toString();
		System.out.println("String of output from result1 is");
		System.out.println(s1);
		
		cr1 = null;
		
		Conrad cr2 = new Conrad(fileModel2);
		cr2.train(fileData);
		double[] weights2 = cr2.getWeights();	
		
		cr2.test(fileData);	
		
		String s2 = cr2.getOutputHandler().toString();
		System.out.println("String of output from result2 is");
		System.out.println(s2);
		
		Assert.a(weights2.length == weights1.length,"length1 is " + weights1.length  + "   and weights2 is " + weights2.length);
		assertArrayEquals(weights2,weights1, tolerance);
		
		assertEquals(s1,s2);
	}	

	public void compareTwoEquivalentModelsTrainTestSameDataWithDifferentStates(String fileModel1, String fileModel2, String fileData, double tolerance ) throws IOException {
		Conrad cr1 = new Conrad(fileModel1);
		cr1.train(fileData);
		double[] weights1 = cr1.getWeights();	
		
		cr1.test(fileData);	
		
		String s1 = cr1.getOutputHandler().toString();
		System.out.println("String of output from result1 is");
		System.out.println(s1);
		
		cr1 = null;
		
		Conrad cr2 = new Conrad(fileModel2);
		cr2.train(fileData);
		double[] weights2 = cr2.getWeights();	
		
		cr2.test(fileData);	
		
		String s2 = cr2.getOutputHandler().toString();
		System.out.println("String of output from result2 is");
		System.out.println(s2);
		
		Assert.a(weights2.length == weights1.length,"length1 is " + weights1.length  + "   and weights2 is " + weights2.length);
		assertArrayEquals(weights2,weights1, tolerance);
	}
	
}
