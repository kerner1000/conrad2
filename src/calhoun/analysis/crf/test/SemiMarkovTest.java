package calhoun.analysis.crf.test;

import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.Conrad;
import calhoun.analysis.crf.CRFInference.InferenceResult;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.util.AbstractTestCase;
import calhoun.util.Assert;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

public class SemiMarkovTest extends AbstractTestCase {
	@SuppressWarnings("unused")
	private static final Log log = LogFactory.getLog(SemiMarkovTest.class);
	
	public void testSemiMarkovBadTraining() throws Exception {
		// Assert that trying to train the semi-Markov fails when training data is used that contains segments longer than the maximum allowed length
		String failureMessage = null;
		try {
			Conrad.main(new String[] {"train", "test/input/zeroOrderLBFGSCachedSemiMarkov.xml", "test/input/zeroOrderTest.txt", "test/working/zeroLBGFSModelCachedSemiMarkov.ser"});
		}
		catch(Exception ex) {
			failureMessage = ex.getMessage();
		}
		assertEquals("Seq #0 Pos 150 Training segment 127 is longer than allowed length 20", failureMessage);
	}

	public void testSemiCRFViterbiCompareWithBaseClass() throws Exception {
		Conrad r = new Conrad("test/input/semiMarkovTestModelNoExplicitLengths.xml");
		r.trainFeatures("test/input/zeroOrderTrivial.txt");
		r.setWeights(new double[] {1,1,1});

		Conrad s = new Conrad("test/input/semiMarkovTestModelNoExplicitLengthsUseBaseClass.xml");
		s.trainFeatures("test/input/zeroOrderTrivial.txt");
		s.setWeights(new double[] {1,1,1});
		
		doViterbiComparison("test/input/zeroOrderTrivial.txt", r, s);
		doViterbiComparison("test/input/zeroOrderTest.txt", r, s);
	}
	
	public void testSemiCRFViterbiCompareWithBaseClassWithLengths() throws Exception {
		Conrad r = new Conrad("test/input/semiMarkovTestModelNoExplicitLengthFeatures.xml");
		r.trainFeatures("test/input/zeroOrderTrivial.txt");
		r.setWeights(new double[] {1,1,1});

		Conrad s = new Conrad("test/input/semiMarkovTestModelNoExplicitLengthsUseBaseClass.xml");
		s.trainFeatures("test/input/zeroOrderTrivial.txt");
		s.setWeights(new double[] {1,1,1});
		
		doViterbiComparison("test/input/zeroOrderTrivial.txt", r, s);
	}
	
	public void testSemiCRFViterbiCompareWithBaseClassWithFeatures() throws Exception {
		Conrad r = new Conrad("test/input/semiMarkovTestModelHalfAndHalf.xml");
		r.trainFeatures("test/input/zeroOrderTrivial.txt");
		r.setWeights(new double[] {1,1,1});

		Conrad s = new Conrad("test/input/semiMarkovTestModelNoExplicitLengthsUseBaseClass.xml");
		//Conrad s = new Conrad("test/input/semiMarkovTestModelNoExplicitLengthFeatures.txt");
		s.trainFeatures("test/input/zeroOrderTrivial.txt");
		s.setWeights(new double[] {1,1,1});
		
		doViterbiComparison("test/input/zeroOrderTrivial.txt", r, s);
	}

	public void testSemiCRFCompareWithBaseClass() throws Exception {
		Conrad s = new Conrad("test/input/semiMarkovTestModelNoExplicitLengthsUseBaseClass.xml");
		s.train("test/input/zeroOrderTrivial.txt");
		
		Conrad r = new Conrad("test/input/semiMarkovTestModelNoExplicitLengths.xml");
		r.train("test/input/zeroOrderTrivial.txt");
		assertEquals(r.getWeights()[0], s.getWeights()[0], 0.0001);

		r = new Conrad("test/input/semiMarkovTestModelNoExplicitLengthFeatures.xml");
		r.train("test/input/zeroOrderTrivial.txt");
		assertEquals(r.getWeights()[0], s.getWeights()[0], 0.0001);

		r = new Conrad("test/input/semiMarkovTestModelHalfAndHalf.xml");
		r.train("test/input/zeroOrderTrivial.txt");
		assertEquals(r.getWeights()[0], s.getWeights()[0], 0.001);
	}
	
	void doViterbiComparison(String file, Conrad a, Conrad b) throws Exception {
		double[] M,N;
		
		List<? extends TrainingSequence<?>> train = a.getInputHandler().readTrainingData(file);
		for(TrainingSequence<?> seq : train) {
			InferenceResult r1 = a.predict(seq.getInputSequence());
			InferenceResult r2 = b.predict(seq.getInputSequence());
			M = r1.bestScores;
			N = r2.bestScores;
			Assert.a(M.length == a.getNumStates());
			
			// Scores should be identical
			for(int r = 0; r<M.length; ++r) {
					assertEquals("State: "+r+" ", N[r], M[r], 0.0001);
			}
			
			// Final paths are not identical due to floating point rounding errors
			//for(int i=0; i < r1.hiddenStates.length; ++i) {
			//	assertEquals("Difference at pos "+i, r1.hiddenStates[i], r2.hiddenStates[i]);
			//}
		}
	}
	
	public void testColt() {
		DenseDoubleMatrix2D f = new DenseDoubleMatrix2D(2,2);
		DenseDoubleMatrix2D g = new DenseDoubleMatrix2D(2,2);
		f.setQuick(1,1,3);
		g.setQuick(1,1,3);
		assertEquals(f,g);
	}
}
