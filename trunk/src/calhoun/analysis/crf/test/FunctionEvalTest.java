package calhoun.analysis.crf.test;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.CRFObjectiveFunctionGradient;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.SemiMarkovSetup;
import calhoun.analysis.crf.io.IntInput;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.analysis.crf.solver.CacheProcessorDeluxe;
import calhoun.analysis.crf.solver.MaximumLikelihoodGradient;
import calhoun.analysis.crf.solver.MaximumLikelihoodSemiMarkovGradient;
import calhoun.analysis.crf.solver.check.AllSparseLengthCacheProcessor;
import calhoun.analysis.crf.solver.check.BasicCRFGradient;
import calhoun.analysis.crf.solver.check.NormalizedCRFGradient;
import calhoun.util.AbstractTestCase;
import calhoun.util.ColtUtil;

public class FunctionEvalTest extends AbstractTestCase {
	private static final Log log = LogFactory.getLog(FunctionEvalTest.class);

	// Test function and gradient calculation using just edge features
	public void testGradientEdgeTrivial() throws Exception {
		doLikelihoodTest(1, IntInput.prepareData("0\n0"), -.6931, 0.000, -.6931, 0.000);
	}

	// / Test function and gradient calculation using just edge features
	public void testGradientEdgeShort() throws Exception {
		double a = Math.log(.5 * 2 / 3.0);
		double b = Math.log(.5 * 2 / 3.0 * 2 / 3.0);
		// grad = 1/3 ln 2 = .231
		doLikelihoodTest(1, IntInput.prepareData("00\n00"), a/2.0, 0.2310/2.0, -.9163/2.0, 0.1386/2.0);
		// grad = 2/3 ln 2 = .462
		doLikelihoodTest(1, IntInput.prepareData("000\n000"), b/3.0, 0.4620/3.0);
	}

	// Test function and gradient calculation using just node features
	public void testGradStateFunc() throws Exception {
		doLikelihoodTest(0, IntInput.prepareData("001111\n001111"), -3.8191/6.0, 0.000/6.0, -4.1115/6.0, -.55451/6.0);
	}

	// Test 2 features
	public void testTwoFeaturesTrivial() throws Exception {
		doLikelihoodTest(2, IntInput.prepareData("00\n00"), -1.9459/2.0, -.990/2.0);
		doLikelihoodTest(2, IntInput.prepareData("00\n00\n00\n00"), -3.8918/4.0, -1.9804/4.0);
	}

	// Test 2 features non-trivial
	public void testTwoFeaturesNonTrivial() throws Exception {
		doLikelihoodTest(2,
				IntInput.prepareData("00001010100100111000\n00001010100100111000\n00001010100100111001\n00001010100100111001\n"), -39.0440/40.0, -11.2550/40.0);
	}

	void doLikelihoodTest(int mmNum, List<? extends TrainingSequence<?>> data, double f, double g) {
		doLikelihoodTest(mmNum, data, f, g, false);
	}

	void doLikelihoodTest(int mmNum, List<? extends TrainingSequence<?>>  data, double f, double g, double h, double i) {
		doLikelihoodTest(mmNum, data, f, g, false);
		doLikelihoodTest(mmNum, data, h, i, true);
	}

	void doLikelihoodTest(int mmNum, List<? extends TrainingSequence<?>>  data, double f, double g, boolean skewedWeights) {
		ModelManager m = new TestFeatureManager(mmNum);
		double[] weights = skewedWeights ? new double[] { 2, 0.5 } : new double[] { 1, 1 };
		double[] grad = new double[2];

		CRFObjectiveFunctionGradient gradFunc;
		double val;
		AllSparseLengthCacheProcessor cacheProcessor;
		short[] max = new short[2];
		
		gradFunc= new BasicCRFGradient();
		gradFunc.setTrainingData(m, data); 
		val = gradFunc.apply(weights, grad);
		log.info("Grad(Basic): " + ColtUtil.format(grad));
		assertEquals(f, val, 0.001);
		assertEquals(g, grad[0], 0.001);

		gradFunc = new NormalizedCRFGradient();
		gradFunc.setTrainingData(m, data); 
		val = gradFunc.apply(weights, grad);
		log.info("Grad(Norm): " + ColtUtil.format(grad));
		assertEquals(f, val, 0.001);
		assertEquals(g, grad[0], 0.001);

		cacheProcessor = new AllSparseLengthCacheProcessor();
		cacheProcessor.setAllPaths(true);
		
		gradFunc = new MaximumLikelihoodGradient();
		cacheProcessor = new AllSparseLengthCacheProcessor();
		cacheProcessor.setAllPaths(true);
		((MaximumLikelihoodGradient) gradFunc).setCacheProcessor(cacheProcessor);
		gradFunc.setTrainingData(m, data); 
		val = gradFunc.apply(weights, grad);
		log.info("Grad(Cache,Valid Paths): " + ColtUtil.format(grad));
		assertEquals(f, val, 0.001);
		assertEquals(g, grad[0], 0.001);

		gradFunc = new MaximumLikelihoodGradient();
		cacheProcessor = new AllSparseLengthCacheProcessor();
		cacheProcessor.setAllPaths(false);
		((MaximumLikelihoodGradient) gradFunc).setCacheProcessor(cacheProcessor);
		gradFunc.setTrainingData(m, data); 
		val = gradFunc.apply(weights, grad);
		log.info("Grad(Cache): " + ColtUtil.format(grad));
		assertEquals(f, val, 0.001);
		assertEquals(g, grad[0], 0.001);

		// Semi CRF no lookback
		Arrays.fill(max, (short) 1);
		gradFunc = new MaximumLikelihoodSemiMarkovGradient();
		cacheProcessor = new AllSparseLengthCacheProcessor();
		cacheProcessor.setAllPaths(false);
		cacheProcessor.setSemiMarkovSetup(new SemiMarkovSetup(max));
		((MaximumLikelihoodSemiMarkovGradient) gradFunc).setCacheProcessor(cacheProcessor);
		gradFunc.setTrainingData(m, data); 

		val = gradFunc.apply(weights, grad);
		log.info("Grad(Semi-cache): " + ColtUtil.format(grad));
		assertEquals(f, val, 0.001);
		assertEquals(g, grad[0], 0.001);

		// Semi CRF lookback without length features
		// We need a lookback of 20 because otherwise we can't have an equivalent test between semi-Markov and regular features here.
		max = new short[2];
		Arrays.fill(max, (short) 20);
		m = new TestFeatureManager(mmNum, 1);
		/*gradFunc = new CachedSemiCRFGradient(max, false);
		gradFunc.setTrainingData(m, data); 
		val = gradFunc.apply(weights, grad);
		log.info("Grad(Length): " + ColtUtil.format(grad));
		assertEquals(f, val, 0.001);
		assertEquals(g, grad[0], 0.001);*/

		gradFunc = new MaximumLikelihoodSemiMarkovGradient();
		CacheProcessorDeluxe cacheProcessor1 = new CacheProcessorDeluxe();
		cacheProcessor1.setAllPaths(false);
		cacheProcessor1.setSemiMarkovSetup(new SemiMarkovSetup(max, true));
		((MaximumLikelihoodSemiMarkovGradient) gradFunc).setCacheProcessor(cacheProcessor1);
		gradFunc.setTrainingData(m, data); 
		val = gradFunc.apply(weights, grad);
		log.info("Grad(Length): " + ColtUtil.format(grad));
		assertEquals(f, val, 0.001);
		assertEquals(g, grad[0], 0.001);

		// Semi CRF lookback with length features
		m = new TestFeatureManager(mmNum, 2);
		Arrays.fill(max, (short) 20);
		gradFunc = new MaximumLikelihoodSemiMarkovGradient();
		cacheProcessor = new AllSparseLengthCacheProcessor();
		cacheProcessor.setAllPaths(false);
		cacheProcessor.setSemiMarkovSetup(new SemiMarkovSetup(max, true));
		//((MaximumLikelihoodSemiMarkovGradient) gradFunc).setBetaLengthFile("betaSM.txt");
		((MaximumLikelihoodSemiMarkovGradient) gradFunc).setCacheProcessor(cacheProcessor);
		gradFunc.setTrainingData(m, data); 
		val = gradFunc.apply(weights, grad);
		log.info("Grad(Length): " + ColtUtil.format(grad));
		assertEquals(f, val, 0.001);
		assertEquals(g, grad[0], 0.001);
	}

}
