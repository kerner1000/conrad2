package calhoun.analysis.crf.test;

import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.Conrad;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.SemiMarkovSetup;
import calhoun.analysis.crf.CacheStrategySpec.CacheStrategy;
import calhoun.analysis.crf.io.IntInput;
import calhoun.analysis.crf.io.StringInput;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.analysis.crf.solver.CacheProcessor;
import calhoun.analysis.crf.solver.CacheProcessorDeluxe;
import calhoun.analysis.crf.solver.MaximumLikelihoodSemiMarkovGradient;
import calhoun.analysis.crf.solver.NoCachingCacheProcessor;
import calhoun.analysis.crf.solver.StandardOptimizer;
import calhoun.analysis.crf.solver.CacheProcessor.FeatureEvaluation;
import calhoun.analysis.crf.solver.CacheProcessor.LengthFeatureEvaluation;
import calhoun.analysis.crf.solver.check.AllSparseLengthCacheProcessor;
import calhoun.util.AbstractTestCase;
import calhoun.util.Assert;

public class CacheProcessorTest extends AbstractTestCase {
	private static final Log log = LogFactory.getLog(CacheProcessorTest.class);

	public void testCPDRejectTrainingDataStatesTooLong() throws Exception {
		checkFailure("test/input/interval13/config/shortIntergenicCPD.xml", "test/input/interval13/data/tooLong.txt");
	}
	
	public void testCPDRejectTrainingDataStatesTooShort() throws Exception {
		checkFailure("test/input/interval13/config/lengthDependentCPDuncommented.xml", "test/input/interval13/data/tooShort.txt");
	}
	
	public void testAllSparseRejectTrainingDataStatesTooShort() throws Exception {
		checkFailure("test/input/interval13/config/lengthDependentAllSparse.xml", "test/input/interval13/data/tooShort.txt");
	}
	
	public void testCPDRejectTrainingDataStatesTooShortStart() throws Exception {
		checkFailure("test/input/interval13/config/lengthDependentCPDuncommented.xml", "test/input/interval13/data/tooShortStart.txt");
	}
	
	public void testCPDRejectTrainingDataStatesTooShortEnd() throws Exception {
		checkFailure("test/input/interval13/config/lengthDependentCPDuncommented.xml", "test/input/interval13/data/tooShortEnd.txt");
	}
	
	public void testCPDRejectTrainingDataStatesViolatesConstraints() throws Exception {
		checkFailure("test/input/interval13/config/lengthDependentCPDuncommented.xml", "test/input/interval13/data/badConstraints.txt");
	}
	
	public void testCPDAllSequencesInvalid() throws Exception {
		checkFailure("test/input/interval13/config/lengthDependentCPDuncommentedDiscard.xml", "test/input/interval13/data/tooShortEnd.txt");
	}
	
	public void testCPDAllDiscardInvalid() throws Exception {
		Conrad conrad = new Conrad("test/input/interval13/config/lengthDependentCPDuncommentedDiscard.xml");
		conrad.train("test/input/interval13/data/oneGoodOneBad.txt");
	}
	
	public void testCPDAllDiscardInvalidLocalScore() throws Exception {
		Conrad conrad = new Conrad("test/input/interval13/config/lengthDependentCPDuncommentedDiscardLocalScore.xml");
		conrad.train("test/input/interval13/data/oneGoodOneBad.txt");
	}
	
	void checkFailure(String configFile, String data) {
		Conrad conrad = new Conrad(configFile);
		boolean fail = false;
		try {
			conrad.trainFeatures(data);
			conrad.trainWeights(conrad.getInputHandler().readTrainingData(data));
		}
		catch(Exception ex) {
			log.warn(ex);
			fail = true;
		}
		assertTrue(fail);
	}
	
	// Test edge features with 1 node, should be no evaluations
	public void testEdgeTrivial() throws Exception {
		int[][] indices = new int[][] { {-1}, {-1}, {-1}, {-1}, {-1}, {-1} };
		float[][] vals = new float[0][0];
		double[] featureSums = new double[] { 0.0 };
		List<? extends TrainingSequence<?>> data = IntInput.prepareData("0\n0"); 
		System.out.println("number of data sequences is " + data.size());
		doTest(1, data, 0, 0, indices, vals, featureSums);
	}

	// Test edge features with 2 and 3 positions.
	public void testEdgeShort() throws Exception {
		int[][] indices = new int[][] { {-1}, {-1}, {0, -1}, {0, -1}, {0, -1}, {0, -1} };
		float[][] vals = new float[][] { {}, {}, {-0.4054651f}, {-1.0986123f}, {-1.0986123f}, {-0.4054651f} };
		double[] featureSums = new double[] { -0.4054651f };
		doTest(1, IntInput.prepareData("00\n00"), 0, 1, indices, vals, featureSums);
	}

	// Test node features
	public void testNode() throws Exception {
		int[][] indices = new int[][] { {0, -1}, {0, -1}, {-1}, {-1}, {-1}, {-1} };
		float[][] vals1 = new float[][] { {-1.0986123f}, {-0.4054651f}, {}, {}, {}, {} };
		double[] featureSums = new double[] { -1.0986123f*2 + -0.4054651f*4 };
		doTest(0, IntInput.prepareData("001111\n001111"), 0, 0, indices, vals1, featureSums); // the DCP feature sums for this are wrong
		float[][] vals2 = new float[][] { {-1.0986123f}, {-0.4054651f}, {}, {}, {}, {} };
		doTest(0, IntInput.prepareData("001111\n001111"), 0, 4, indices, vals2, featureSums);
	}

	// Test 2 features
	public void testTwoFeaturesTrivial() throws Exception {
		int[][] indices1 = new int[][] { {0, -1}, {0, -1}, {-1}, {-1}, {-1}, {-1} };
		float[][] vals1 = new float[][] { {-1.0986123f}, {-0.4054651f}, {}, {}, {}, {} };
		double[] featureSums = new double[] { -1.0986123f*2, -0.4054651f };
		doTest(2, IntInput.prepareData("00\n00"), 0, 0, indices1, vals1, featureSums);
		int[][] indices2 = new int[][] { {0, -1}, {0, -1}, {1, -1}, {1, -1}, {1, -1}, {1, -1} };
		float[][] vals2 = new float[][] { {-1.0986123f}, {-0.4054651f}, {-0.4054651f}, {-1.0986123f}, {-1.0986123f}, {-0.4054651f} };
		featureSums = new double[] { -1.0986123f*4, -0.4054651f*2 };
		doTest(2, IntInput.prepareData("00\n00\n00\n00"), 1, 1, indices2, vals2, featureSums);
	}

	// Test 2 features non-trivial
	public void testTwoFeaturesNonTrivial() throws Exception {
		int[][] indices = new int[][] { {0, -1}, {0, -1}, {1, -1}, {1, -1}, {1, -1}, {1, -1} };
		float[][] vals = new float[][] { {-1.0986123f}, {-0.4054651f}, {-0.4054651f}, {-1.0986123f}, {-1.0986123f}, {-0.4054651f} };
		double[] featureSums = new double[] { -33.54728, -29.963765};
		doTest(2, IntInput.prepareData("00001010100100111000\n00001010100100111000\n00001010100100111001\n00001010100100111001\n"), 1, 4, indices, vals, featureSums);
	}

	public void testLengthCacheDummy() throws Exception {
		int[][] lookbacks = new int[][] { {0, 1, -1}, {0, 1, -1} };
		int[][] nodeIndices = new int[][] { {0, 0, 0, -1}, {1, 0, 0, -1} };
		float[][] nodeValues = new float[0][4];
		ModelManager m = new TestFeatureManager(1);
		doLengthTest(m, IntInput.prepareData("00\n00"), 0, 1, 2, lookbacks, nodeIndices, nodeValues);

		m = new TestFeatureManager(2);
		lookbacks = new int[][] { {0, 1, 2, 3, -1}, {0, 1, 2, 3, -1} };
		nodeIndices = new int[][] { {0, 0, 0, -1}, {0, 1, 0, -1}, {0, 2, 0, -1}, {0, 3, 0, -1}, {1, 0, 0, -1} };
		doLengthTest(m, IntInput.prepareData("0000\n0000"), 0, 3, 2, lookbacks, nodeIndices, nodeValues);
	}

	public void testLengthCache() throws Exception {
		int[][] lookbacks = new int[][] { {0, 1, 2, 3, -1} };
		int[][] nodeIndices = new int[0][4];
		float[][] nodeValues = new float[][] { {0, 0, 0, -0.0f}, {0, 1, 0, -0.11157f} , {0, 2, 0, -0.11157f*2}, {0, 3, 0, -0.11157f*3}};
		Conrad c = new Conrad("test/input/semiMarkovTestModelHalfAndHalf.xml");
		ModelManager m = c.getModel();
		List<? extends TrainingSequence<?>> data = StringInput.prepareData("00110\nATGCA");
		c.trainFeatures(data);
		CacheProcessor cp = ((MaximumLikelihoodSemiMarkovGradient) ((StandardOptimizer)c.getOptimizer()).getObjectiveFunction()).getCacheProcessor();
		cp.setTrainingData(m, data);
		cp.evaluateSegmentsEndingAt(0, 3);
		LengthFeatureEvaluation[][] lenEvals = cp.getLengthFeatureEvaluations();
		checkLengthEvals(lenEvals, 1, lookbacks, nodeIndices, nodeValues);
	}
	
	void doLengthTest(ModelManager m, List<? extends TrainingSequence<?>> data, int seq, int pos, int nStates, int[][] lookback, int[][] nodeIndices, float[][] nodeValues) {
		AllSparseLengthCacheProcessor cp = new AllSparseLengthCacheProcessor();
		SemiMarkovSetup setup = new SemiMarkovSetup(new short[] {4, 4});
		setup.setIgnoreSemiMarkovSelfTransitions(true);
		cp.setSemiMarkovSetup(setup);
		cp.setTrainingData(m, data);
		cp.evaluateSegmentsEndingAt(seq, pos);
		LengthFeatureEvaluation[][] lenEvals = cp.getLengthFeatureEvaluations();
		checkLengthEvals(lenEvals, nStates, lookback, nodeIndices, nodeValues);
	}

	void checkLengthEvals(LengthFeatureEvaluation[][] lenEvals, int nStates, int[][] lookback, int[][] nodeIndices, float[][] nodeValues) {
		assertEquals(nStates, lenEvals.length);
		for(int i=0; i<lookback.length; ++i) {
			for(int j=0; j<lookback[i].length; ++j) {
				assertEquals(lookback[i][j], lenEvals[i][j].lookback);
			}
		}
		for(int[] entry : nodeIndices) {
			assertEquals(entry[3], lenEvals[entry[0]][entry[1]].nodeEval.index[entry[2]]);
		}
		
		for(float[] entry : nodeValues) {
			assertEquals(entry[3], lenEvals[(int)entry[0]][(int)entry[1]].nodeEval.value[(int)entry[2]], 0.0001);
		}
	}
	
	void doTest(int mmNum, List<? extends TrainingSequence<?>> data, int seq, int pos, int[][] indices, float[][] vals, double[] featureSums) {
		ModelManager m = new TestFeatureManager(mmNum);
		
		AllSparseLengthCacheProcessor cp = new AllSparseLengthCacheProcessor();
		testOneCacheProcessor(cp,m,data,seq,pos,indices,vals,featureSums);
		
		NoCachingCacheProcessor ncp = new NoCachingCacheProcessor();
		testOneCacheProcessor(ncp,m,data,seq,pos,indices,vals,featureSums);
	
		CacheProcessorDeluxe dcp = new CacheProcessorDeluxe();
		//dcp.setSemiMarkovSetup(new SemiMarkovSetup(new short[]{1,1}, new short[]{50,50},true));
		testOneCacheProcessor(dcp,m,data,seq,pos,indices,vals,featureSums);

		CacheProcessorDeluxe dcp2 = new CacheProcessorDeluxe(CacheStrategy.CONSTANT);
		//dcp2.setSemiMarkovSetup(new SemiMarkovSetup(new short[]{1,1}, new short[]{50,50},true));
		testOneCacheProcessor(dcp2,m,data,seq,pos,indices,vals,featureSums);
		
		CacheProcessorDeluxe dcp3 = new CacheProcessorDeluxe(CacheStrategy.DENSE);
		//dcp3.setSemiMarkovSetup(new SemiMarkovSetup(new short[]{1,1}, new short[]{50,50},true));
		testOneCacheProcessor(dcp3,m,data,seq,pos,indices,vals,featureSums);
		
		CacheProcessorDeluxe dcp4 = new CacheProcessorDeluxe(CacheStrategy.SPARSE);
		//dcp4.setSemiMarkovSetup(new SemiMarkovSetup(new short[]{1,1}, new short[]{50,50},true));
		testOneCacheProcessor(dcp4,m,data,seq,pos,indices,vals,featureSums);		
	}

	void testOneCacheProcessor(CacheProcessor dcp, ModelManager m, List<? extends TrainingSequence<?>> data,  int seq, int pos, int[][] indices, float[][] vals, double[] featureSums) {
		//dcp.setAllPaths(false);
		dcp.setTrainingData(m, data);
		dcp.evaluatePosition(seq, pos);
		if (pos > 0) {
			assertEvalEquals(dcp.getFeatureEvaluations(), indices, vals);
		} else {
			assertNonedgeEvalEquals(m,dcp.getFeatureEvaluations(), indices, vals);		
		}
		assertArrayEquals(featureSums, dcp.getFeatureSums(), 0.00001);
	}
	
	
	private void assertNonedgeEvalEquals(ModelManager m, FeatureEvaluation[] evals, int[][] indices, float[][] vals) {
		log.warn(evals);
		
		Assert.a(indices.length >= m.getNumStates());
		
		for(int i=0; i<m.getNumStates(); ++i) {
			for(int j=0; j<indices[i].length; ++j) {
				assertEquals(indices[i][j], evals[i].index[j]);
			}
		}
		for(int i=0; i<vals.length; ++i) {
			for(int j=0; j<vals[i].length; ++j) {
				assertEquals(vals[i][j], evals[i].value[j], .00001);
			}
		}
	}
	
	private void assertEvalEquals(FeatureEvaluation[] evals, int[][] indices, float[][] vals) {
		log.warn(evals);
		for(int i=0; i<indices.length; ++i) {
			for(int j=0; j<indices[i].length; ++j) {
				assertEquals(indices[i][j], evals[i].index[j]);
			}
		}
		for(int i=0; i<vals.length; ++i) {
			for(int j=0; j<vals[i].length; ++j) {
				assertEquals(vals[i][j], evals[i].value[j], .00001);
			}
		}
	}
}
