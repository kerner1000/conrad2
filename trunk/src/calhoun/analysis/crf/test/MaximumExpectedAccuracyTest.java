package calhoun.analysis.crf.test;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.Conrad;
import calhoun.analysis.crf.LocalPathSimilarityScore;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.SemiMarkovSetup;
import calhoun.analysis.crf.io.IntInput;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.analysis.crf.scoring.SimScoreMaxStateAgreement;
import calhoun.analysis.crf.solver.CacheProcessorDeluxe;
import calhoun.analysis.crf.solver.MaximumExpectedAccuracySemiMarkovGradient;
import calhoun.analysis.crf.solver.StandardOptimizer;
import calhoun.analysis.crf.solver.check.CachedAOFGradient;
import calhoun.util.AbstractTestCase;
import calhoun.util.ColtUtil;

public class MaximumExpectedAccuracyTest extends AbstractTestCase {
	private static final Log log = LogFactory.getLog(MaximumExpectedAccuracyTest.class);

	public void testLittle() throws Exception {
		
		String config         = "test/input/aofreg_test/delta_conservation_aofreg.xml";
		String input          = "test/input/aofreg_test/testSeq.txt";

		Conrad conrad = new Conrad(config);
		conrad.train(input);
		double[] regularWeights = conrad.getWeights();
		conrad.test(input);

		config         = "test/input/aofreg_test/delta_conservation_aofreg_semi.xml";
		conrad = new Conrad(config);
		conrad.train(input);
		double[] semiMarkovWeights = conrad.getWeights();
		conrad.test(input);
		assertArrayEquals(regularWeights, semiMarkovWeights, 0.001);
	}

	public void testGradEvals() throws Exception {
		doFuncEvalsSkewed(0, IntInput.prepareData("0\n0"));
		doFuncEvalsSkewed(0, IntInput.prepareData("00\n00"));
		doFuncEvalsSkewed(0, IntInput.prepareData("000\n000"));
		doFuncEvalsSkewed(1, IntInput.prepareData("1\n1"));
		doFuncEvalsSkewed(3, IntInput.prepareData("11\n11"));
		doFuncEvalsSkewed(1, IntInput.prepareData("10\n10"));
		doFuncEvalsSkewed(1, IntInput.prepareData("111\n111"));
		doFuncEvalsSkewed(3, IntInput.prepareData("101\n101"));
		doFuncEvalsSkewed(3, IntInput.prepareData("011\n011"));
		doFuncEvalsSkewed(1, IntInput.prepareData("1111\n1111"));
		doFuncEvalsSkewed(1, IntInput.prepareData("001111\n001111"));
		doFuncEvalsSkewed(2, IntInput.prepareData("00\n00"));
		doFuncEvalsSkewed(2, IntInput.prepareData("00\n00\n00\n00"));
		doFuncEvalsSkewed(1, IntInput.prepareData("00\n00\n0\n0\n"));
		doFuncEvalsSkewed(1, IntInput.prepareData("0\n0\n1111\n1111\n"));
		doFuncEvalsSkewed(2, IntInput.prepareData("00001010100100111000\n00001010100100111000\n00001010100100111001\n00001010100100111001\n"));
	}

	void doFuncEvalsSkewed(int mmNum, List<? extends TrainingSequence<?>> data) throws Exception {
		//doFuncEvals(mmNum, false, false, data);
		//doFuncEvals(mmNum, true, false, data);
		
		if(mmNum == 1 || mmNum == 3) {
			//doFuncEvals(mmNum, false, true, data);
			doFuncEvals(mmNum, true, true, data);
		}
	}
	
	void doFuncEvals(int mmNum, boolean skewedWeights, boolean fm3, List<? extends TrainingSequence<?>> data) throws Exception {
		double[] weights = skewedWeights ? new double[] { 2, 0.5, 1 } : new double[] { 1, 1, 1 };

		MaximumExpectedAccuracySemiMarkovGradient gradFunc;
		CacheProcessorDeluxe cacheProcessor;
		short[] max = new short[2];
		Arrays.fill(max, (short) 1);

		ModelManager m = fm3 ? new TestFeatureManager3(mmNum, false) : new TestFeatureManager2(mmNum, false);
		m.train(0, m, data);
		
		// Semi CRF no lookback
		gradFunc = new MaximumExpectedAccuracySemiMarkovGradient();
		cacheProcessor = new CacheProcessorDeluxe();
		cacheProcessor.setAllPaths(false);
		cacheProcessor.setSemiMarkovSetup(new SemiMarkovSetup(max, true));
		((MaximumExpectedAccuracySemiMarkovGradient) gradFunc).setCacheProcessor(cacheProcessor);
		gradFunc.setTrainingData(m, data); 
		double[] grad1 = new double[m.getNumFeatures()];
		double val1 = gradFunc.apply(weights, grad1);
		log.info("Grad1: " + ColtUtil.format(grad1));

		// Semi CRF lookback
		Arrays.fill(max, (short) 20);
		m = fm3 ? new TestFeatureManager3(mmNum, false) : new TestFeatureManager2(mmNum, false);
		m.train(0, m, data);
		gradFunc = new MaximumExpectedAccuracySemiMarkovGradient();
		cacheProcessor = new CacheProcessorDeluxe();
		cacheProcessor.setAllPaths(false);
		cacheProcessor.setSemiMarkovSetup(new SemiMarkovSetup(max, true));
		((MaximumExpectedAccuracySemiMarkovGradient) gradFunc).setCacheProcessor(cacheProcessor);
		//((MaximumExpectedAccuracySemiMarkovGradient) gradFunc).setBetaLengthFile("modeMarginal.txt");
		//((MaximumExpectedAccuracySemiMarkovGradient) gradFunc).setNodeMarginalFile("nodeMarginal.txt");
		//((MaximumExpectedAccuracySemiMarkovGradient) gradFunc).setExpectedProductFile("expectedProduct.txt");
		//((MaximumExpectedAccuracySemiMarkovGradient) gradFunc).setScoreAlphaFile("scoreAlpha.txt");
		gradFunc.setTrainingData(m, data); 
		double[] grad2 = new double[m.getNumFeatures()];
		double val2 = gradFunc.apply(weights, grad2);
		log.info("Grad2: " + ColtUtil.format(grad2));
		
		assertEquals(val1, val2, 0.001);
		assertArrayEquals(grad1, grad2, 0.001);

		// Semi CRF lookback with length features
		m = fm3 ? new TestFeatureManager3(mmNum, true) : new TestFeatureManager2(mmNum, true);
		m.train(0, m, data);
		gradFunc = new MaximumExpectedAccuracySemiMarkovGradient();
		cacheProcessor = new CacheProcessorDeluxe();
		cacheProcessor.setAllPaths(false);
		cacheProcessor.setSemiMarkovSetup(new SemiMarkovSetup(max, true));
		((MaximumExpectedAccuracySemiMarkovGradient) gradFunc).setCacheProcessor(cacheProcessor);
		//((MaximumExpectedAccuracySemiMarkovGradient) gradFunc).setExpectedProductFile("expectedProductLen.txt");
		gradFunc.setTrainingData(m, data); 
		double[] grad3 = new double[m.getNumFeatures()];
		double val3 = gradFunc.apply(weights, grad3);
		log.info("Grad3: " + ColtUtil.format(grad3));
		
		assertEquals(val1, val3, 0.001);
		assertArrayEquals(grad1, grad3, 0.001);
	}
	
	public void testGeneCallerLocalScore() throws Exception {
		Conrad nodeOnly = new Conrad("test/input/geneCallerLocal/baseline_aof.xml");
		Conrad length = new Conrad("test/input/geneCallerLocal/baseline_aof_length.xml");
		nodeOnly.train("test/input/geneCallerLocal");
		length.train("test/input/geneCallerLocal");
		assertArrayEquals(nodeOnly.getWeights(), length.getWeights(), 0.0001);
	}
	
	public void testAlternateObjectiveFunction() throws Exception {
		LocalPathSimilarityScore s = new SimScoreMaxStateAgreement();
		// s measures the number of agreeing nucleotides between a hidden path
		// and the actual hidden path, not counting position zero.
		
		// This is a two feature example I worked out by hand
		ModelManager m2 = new TestFeatureManager(2);
		List<? extends TrainingSequence<?>> data = IntInput.prepareData("00\n00");		
		doAlternateObjectiveFunctionTest( m2, data, s, 0.2857142/2.0,-0.1838962/2.0);
		
		// below should be exactly double what is above (the
		// observations are not used by this feature manager)
		List<? extends TrainingSequence<?>> data2 = IntInput.prepareData("00\n11\n00\n00");
		doAlternateObjectiveFunctionTest( m2, data2, s,0.57142853/4.0,-0.367792/4.0);		
		
		// the following should grow linearly with length-1
		// The manager has just one node feature.  At weight 1.0, the probability of a zero,
		// and hence the expected value of S, is 1/3.
		ModelManager m0 = new TestFeatureManager(0);		// this has just node features
		List<? extends TrainingSequence<?>> data3 = IntInput.prepareData("00\n00");
		double t3 = -0.1540327;    		// I didn't check the gradient by hand but copied it from output
		doAlternateObjectiveFunctionTest( m0, data3, s,0.33333333/2.0,t3/2.0);			

		// since this grows linearly, it should be exactly three times the above
		List<? extends TrainingSequence<?>> data4 = IntInput.prepareData("0000\n0000");
		doAlternateObjectiveFunctionTest( m0, data4, s, 1.0/4.0,3*t3/4.0);	

		// Expected value of S should be 0.5; gradient should be zero
		ModelManager m1 = new TestFeatureManager(1);		// this has just edge features
		doAlternateObjectiveFunctionTest( m1, data3, s, 0.5/2.0,0.0/2.0);	

		// Expected value of S should be 1.5; gradient should be 0
		doAlternateObjectiveFunctionTest( m1, data4, s,1.5/4.0,0.0/4.0);
	}

	// Similar logic to above bug on a slightly bigger example.
	// If LBGFS can converge given my gradient/value function, then my function is
	// probably self-consistent.
	public void testSemiRealExample() throws Exception {
		String input          = "test/input/cryptoAOFUnittest/Tiny_1_Train_Test.txt";
		String config1        = "test/input/cryptoAOFUnittest/delta_aof_Model.xml";
		String config2         = "test/input/cryptoAOFUnittest/delta_aof_Model_semi_nolen.xml";
		String config3         = "test/input/cryptoAOFUnittest/delta_aof_Model_semi.xml";
		doSameWeightsTest(config1, config2, config3, input);
	}

	void doSameWeightsTest(String config1, String config2, String config3, String input) throws Exception {
		double[] semiMarkovWeights;
		Conrad conrad;
		conrad = new Conrad(config1);
		conrad.train(input);
		double[] regularWeights = conrad.getWeights();
		conrad.test(input);

		conrad = new Conrad(config2);
		conrad.train(input);
		double[] noLenWeights = conrad.getWeights();
		conrad.test(input);
		assertArrayEquals(regularWeights, noLenWeights, 0.001);

		if(config3 != null) {
			conrad = new Conrad(config3);
			conrad.train(input);
			semiMarkovWeights = conrad.getWeights();
			conrad.test(input);
			assertArrayEquals(noLenWeights, semiMarkovWeights, 0.001);
		}
	}
	
	// The idea of the test below is that of requiring LBGFS to converge to a point with zero gradient.
	// Doesn't this require that the gradient equal zero at a point where the function is maximized.
	// Isn't this very unlikely to succeed if either the function value or its
	// gradient were computed incorrcetly?  I think actually this is a strong test;
	// the gradient defines a direction for a line search, and the function value determines
	// where on that line search to stop next.  It's hard to image this process
	// converging in 2 or more dimensions unless the function and gradient are matched.
	public void testCachedAOFGradient() throws Exception {
		List<? extends TrainingSequence<?>> data = IntInput.prepareData("001111\n001111\n001111\n001111\n001111\n001111\n001111\n001111\n");
		ModelManager m = new TestFeatureManager(2);

		StandardOptimizer opt = new StandardOptimizer();
		opt.setStarts(new double[] {0.1, 0.2});
		opt.setRequireConvergence(true);
		opt.setEpsForConvergence(0.0000005);
		opt.setObjectiveFunction(new CachedAOFGradient());
		opt.optimize(m, data);

		MaximumExpectedAccuracySemiMarkovGradient semiAof = new MaximumExpectedAccuracySemiMarkovGradient();
		semiAof.setCacheProcessor(new CacheProcessorDeluxe());
		opt.setObjectiveFunction(semiAof);
		opt.optimize(m, data);
	}

	void doAlternateObjectiveFunctionTest(ModelManager m, List<? extends TrainingSequence<?>> data, LocalPathSimilarityScore s, double hand_val, double hand_grad0 ) {

		double[] weights = new double[m.getNumFeatures()];
		Arrays.fill(weights,1.0);

		CachedAOFGradient gradFunc = new CachedAOFGradient();
		gradFunc.setAllPaths(true);
		gradFunc.setScoreAlphaFile("scoreAlphaOld.txt");
		gradFunc.setExpectedProductFile("expectedProductOld.txt");
		gradFunc.setLocalPathSimilarityScore(s);
		gradFunc.setTrainingData(m, data);

		double val;
		double[] grad = new double[m.getNumFeatures()];
		
		val = gradFunc.apply(weights, grad);
		log.info("  val = " + val + "   grad[0] = " + grad[0]) ;
		log.info("Grad(Cache,Valid Paths): " + ColtUtil.format(grad));
		assertEquals(hand_val, val, 0.001);
		assertEquals(hand_grad0, grad[0], 0.001);		

		MaximumExpectedAccuracySemiMarkovGradient semiAof = new MaximumExpectedAccuracySemiMarkovGradient();
		semiAof.setCacheProcessor(new CacheProcessorDeluxe());
		semiAof.setTrainingData(m, data);
		semiAof.setMarginalsFile("marginals.txt");
		semiAof.setScoreAlphaFile("scoreAlpha.txt");
		semiAof.setExpectedProductFile("expectedProduct.txt");
		val = semiAof.apply(weights, grad);
		log.info("  val = " + val + "   grad[0] = " + grad[0]) ;
		log.info("Grad(Cache,Valid Paths): " + ColtUtil.format(grad));
		assertEquals(hand_val, val, 0.001);
		assertEquals(hand_grad0, grad[0], 0.001);		
	}
}
