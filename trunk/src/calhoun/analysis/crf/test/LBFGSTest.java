package calhoun.analysis.crf.test;

import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.Conrad;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.io.IntInput;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.analysis.crf.solver.MaximumLikelihoodGradient;
import calhoun.analysis.crf.solver.StandardOptimizer;
import calhoun.util.AbstractTestCase;

public class LBFGSTest extends AbstractTestCase {
	@SuppressWarnings("unused")
	private static final Log log = LogFactory.getLog(LBFGSTest.class);

	/// Simple model, 2 states, 1 feature that returns log 1/3 and log 2/3 for every position.  With weights 1, normalization should be 1.
	public void testLBFGS() throws Exception {
		List<? extends TrainingSequence<?>> data = IntInput.prepareData("001111\n001111\n001111\n001111\n001111\n001111\n001111\n001111\n");
		ModelManager m = new TestFeatureManager(1);

		StandardOptimizer opt = new StandardOptimizer();
		opt.setStarts(new double[] {0.1});
		opt.setEpsForConvergence(0.001);
		opt.setObjectiveFunction(new MaximumLikelihoodGradient());
		double[] weights = opt.optimize(m, data);
		assertEquals(2.0, weights[0], 0.01);
	}

	public void testMainLBGFS() throws Exception {
		// Create and train the CRF
		Conrad.main(new String[] {"train", "test/input/zeroOrderLBFGS.xml", "test/input/zeroOrderTest.txt", "test/working/zeroLBGFSModel.ser"});
		Conrad.main(new String[] {"test", "test/working/zeroLBGFSModel.ser", "test/input/zeroOrderTest.txt", "test/working/zeroOrderMainLBGFSPredicted.txt"});
		Conrad.main(new String[] {"test", "test/working/zeroLBGFSModel.ser", "test/input/zeroOrderTest2.txt", "test/working/zeroOrder2MainLBGFSPredicted.txt"});
		assertFilesMatch("test/output/zeroOrder2MainLBGFSPredicted.txt", "test/working/zeroOrder2MainLBGFSPredicted.txt");
	}

	public void testMainLBGFSCached() throws Exception {
		// Create and train the CRF
		Conrad.main(new String[] {"train", "test/input/zeroOrderLBFGSCached.xml", "test/input/zeroOrderTest.txt", "test/working/zeroLBGFSModelCached.ser"});
		Conrad.main(new String[] {"test", "test/working/zeroLBGFSModelCached.ser", "test/input/zeroOrderTest.txt", "test/working/zeroOrderMainLBGFSCachedPredicted.txt"});
		Conrad.main(new String[] {"test", "test/working/zeroLBGFSModelCached.ser", "test/input/zeroOrderTest2.txt", "test/working/zeroOrder2MainLBGFSCachedPredicted.txt"});
		assertFilesMatch("test/output/zeroOrder2MainLBGFSPredicted.txt", "test/working/zeroOrder2MainLBGFSCachedPredicted.txt");
	}

	public void testMainLBGFSCachedSemiMarkov() throws Exception {
		// Create and train the CRF
		Conrad.main(new String[] {"train", "test/input/zeroOrderLBFGSCachedSemiMarkov.xml", "test/input/zeroOrderTestShortMax.txt", "test/working/zeroLBGFSModelCachedSemiMarkov.ser"});
		Conrad.main(new String[] {"test", "test/working/zeroLBGFSModelCachedSemiMarkov.ser", "test/input/zeroOrderTest.txt", "test/working/zeroOrderMainLBGFSCachedSemiMarkovPredicted.txt"});
		Conrad.main(new String[] {"test", "test/working/zeroLBGFSModelCachedSemiMarkov.ser", "test/input/zeroOrderTest2.txt", "test/working/zeroOrder2MainLBGFSCachedSemiMarkovPredicted.txt"});
		assertFilesMatch("test/output/zeroOrder2MainLBGFSCachedSemiMarkovPredicted.txt", "test/working/zeroOrder2MainLBGFSCachedSemiMarkovPredicted.txt");
	}
}
