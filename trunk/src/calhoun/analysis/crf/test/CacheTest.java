package calhoun.analysis.crf.test;

import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.Conrad;
import calhoun.analysis.crf.io.InputHandlerInterleaved;
import calhoun.analysis.crf.io.StringInput;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.analysis.crf.solver.check.FeatureCache;
import calhoun.util.AbstractTestCase;

public class CacheTest extends AbstractTestCase {
	private static final Log log = LogFactory.getLog(CacheTest.class);
	boolean debug = log.isDebugEnabled();
	
	public void testCache() throws Exception {
		doCache("test/input/zeroOrderModel.xml", "test/input/zeroOrderTrivial.txt", 0, 4);
		doCache("test/input/constraintModel.xml", "test/input/trivialGenes.txt", 0, 29);
		doCache("test/input/semiMarkovTestModelNoExplicitLengths.xml", "test/input/zeroOrderTrivial.txt", 30, 4);
	}
	
	public void testCacheHarder() throws Exception {
		List<? extends TrainingSequence<?>> data = new InputHandlerInterleaved(new StringInput()).readTrainingData("test/input/zeroOrderTest.txt");
		Conrad crf = ZeroOrderManager.getCRF();
		crf.trainFeatures(data);
		
		FeatureCache f = new FeatureCache(crf.getModel(), data);
		assertEquals(798, f.cachedFeatures);
		assertEquals(4, f.numConstantFeatures);
		assertEquals(802, f.totalFeatures);
	}
	
	public void doCache(String model, String data, int cachedFeatures, int constantFeatures) throws Exception {
		Conrad r = new Conrad(model);
		//CRFRunner r = new CRFRunner(model);
		r.trainFeatures(data);
		List<? extends TrainingSequence<?>> training = r.getInputHandler().readTrainingData(data);
		FeatureCache f = new FeatureCache(r.getModel(), training);
		assertEquals(cachedFeatures, f.cachedFeatures);
		assertEquals(constantFeatures, f.numConstantFeatures);
		assertEquals(cachedFeatures + constantFeatures, f.totalFeatures);
	}
}
