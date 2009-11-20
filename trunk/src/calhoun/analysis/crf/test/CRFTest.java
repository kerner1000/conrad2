package calhoun.analysis.crf.test;

import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.Conrad;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.analysis.crf.io.CompositeInput.LegacyInputHandler;
import calhoun.analysis.crf.solver.StandardOptimizer;
import calhoun.util.AbstractTestCase;

public class CRFTest extends AbstractTestCase {
	private static final Log log = LogFactory.getLog(CRFTest.class);

	public void testFullCycle() throws Exception {
		// Create and train the CRF
		Conrad crf = ZeroOrderManager.getCRF();
		List<? extends TrainingSequence<?>> data = crf.getInputHandler().readTrainingData("test/input/zeroOrderTest.txt");
		assertEquals(data.size(),3);

		crf.train(data);
		// Write the trained model out to a file and read it back in again.

		log.info(crf.printWeights());
		crf.test(data, "test/working/zeroOrderPredicted.txt");
		assertFilesMatch("test/output/zeroOrderPredicted.txt", "test/working/zeroOrderPredicted.txt");
		//log.info(results);
	}
	
	public void testSpringConfig() throws Exception {
		Conrad c = new Conrad("test/input/configTest.xml");
		ModelManager cm = c.getModel();
		assertEquals(7, cm.getNumStates());
		assertEquals("exon1", cm.getStateName(1));
		LegacyInputHandler seq = (LegacyInputHandler) c.getInputHandler();
		assertNotNull(seq);
		StandardOptimizer opt = (StandardOptimizer) c.getOptimizer();
		assertEquals(1, opt.getMaxIters());
	}
	
	public void testToolkitMainSeparateTraining() throws Exception {
		// Create and train the CRF
		Conrad.main(new String[] {"trainFeatures", "test/input/zeroOrderLBFGS.xml", "test/input/zeroOrderTest.txt", "test/working/zeroModelSep.ser"});
		Conrad.main(new String[] {"trainWeights", "test/working/zeroModelSep.ser", "test/input/zeroOrderTest.txt", "test/working/zeroModelSep2.ser"});
		Conrad.main(new String[] {"test", "test/working/zeroModelSep2.ser", "test/input/zeroOrderTest.txt", "test/working/zeroOrderSepMainPredicted.txt"});
	}

}
