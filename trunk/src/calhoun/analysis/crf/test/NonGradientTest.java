package calhoun.analysis.crf.test;

import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.CRFObjectiveFunctionGradient;
import calhoun.analysis.crf.Conrad;
import calhoun.analysis.crf.io.InputHandlerInterleaved;
import calhoun.analysis.crf.io.IntInput;
import calhoun.analysis.crf.io.OutputHandlerGeneCallStats;
import calhoun.analysis.crf.io.StringInput;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.analysis.crf.solver.MaximumLikelihoodGradient;
import calhoun.analysis.crf.solver.SimplexOptimizer;
import calhoun.analysis.crf.solver.StandardOptimizer;
import calhoun.analysis.crf.solver.Viterbi;
import calhoun.util.AbstractTestCase;

public class NonGradientTest extends AbstractTestCase {
	static final Log log = LogFactory.getLog(NonGradientTest.class);
	boolean debug = log.isDebugEnabled();
	
	public void testNoGradient() throws Exception {
		SimplexOptimizer opt = new SimplexOptimizer();
		opt.setStepSize(0.000001);
		List<? extends TrainingSequence<?>> data = IntInput.prepareData("001111\n001111");
		opt.setObjectiveFunction(new MaximumLikelihoodGradient());
		double[] weights = opt.optimize(new TestFeatureManager(0), data);
		assertEquals(1.0, weights[0], 0.00001);
	}

	public void testConvergence() throws Exception {
		SimplexOptimizer opt = new SimplexOptimizer();
		opt.setStepSize(0.000001);
		opt.setMaxIters(2);
		opt.setRequireConvergence(true);
		opt.setObjectiveFunction(new MaximumLikelihoodGradient());
		List<? extends TrainingSequence<?>> data = IntInput.prepareData("001111\n001111");
		boolean ex = false;
		try {
			opt.optimize(new TestFeatureManager(0), data);
		}
		catch(Exception e) {
			ex = true;
		}
		assertTrue(ex);
	}

	public void testGradientEdge() throws Exception {
		List<? extends TrainingSequence<?>> data = IntInput.prepareData("00110\n00110");
		CRFObjectiveFunctionGradient obj = new MaximumLikelihoodGradient();

		SimplexOptimizer opt = new SimplexOptimizer();
		opt.setObjectiveFunction(obj);
		double[] weightsSimplex = opt.optimize(new TestFeatureManager(1), data);

		StandardOptimizer stdOpt = new StandardOptimizer(); 
		stdOpt.setObjectiveFunction(obj);
		double[] weightsStandard = stdOpt.optimize(new TestFeatureManager(1), data);
		assertArrayEquals(weightsStandard, weightsSimplex, 1e-5);
	}

	public void testViterbi() throws Exception {
		double[] weights = new double[] { -9.639248, 9.659248, -9.982040e+00, 1.000204e+01, -1.534077e+01, 2.854132e+01, -1.301635e+01, -3.040209e+00, -3.383001e+00 , 1.947956e+01 };
		Conrad crf = new Conrad();
		crf.setModel(new ZeroOrderManager());
		crf.setInputHandler(new InputHandlerInterleaved(new StringInput()));
		OutputHandlerGeneCallStats stats = new OutputHandlerGeneCallStats(crf.getModel(), crf.getInputHandler());
		stats.setWriteTrainingData(true);
		crf.setOutputHandler(stats);
		crf.setInference(new Viterbi());

		crf.trainFeatures("test/input/zeroOrderTest.txt");
		crf.setWeights(weights);
		crf.test("test/input/zeroOrderTest.txt", "test/working/zeroOrderPredicted.txt");
		assertFilesMatch("test/output/zeroOrderPredicted.txt", "test/working/zeroOrderPredicted.txt");
		//log.info(results);
	}
}
