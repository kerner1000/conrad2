package calhoun.analysis.crf.test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.CRFObjectiveFunctionGradient;
import calhoun.analysis.crf.Conrad;
import calhoun.analysis.crf.io.StringInput;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.analysis.crf.solver.StandardOptimizer;
import calhoun.util.AbstractTestCase;

/** Tests that CRF is working with valid probabilities - the sum of all possible labelings is 1.
 *  
 * Test that the code to walk through only the valid paths works correctly.
 * Uses a two state model that disallows transitions to self.  010101... or 101010... are the only allowed paths. */
public class ValidPathsTest extends AbstractTestCase {
	private static final Log log = LogFactory.getLog(ValidPathsTest.class);

	static List<? extends TrainingSequence<?>> train1 = StringInput.prepareData("000001\nATAGAG\n101010\nAGAGAG\n");
	static String[] test1 = new String[] {
			"1\nA\n0\nT\n",
			"10\nAT\n01\nAA\n00\nTA\n11\nAA\n"
	};

	static String[] test2 = new String[] {
		"1\nA\n0\nT\n",
		"10\nAT\n01\nAA\n"
	};

	/** Tests all 2-state combinations with a restricted model and onlyValidTransitions = false */
	public void testSumPathsAllPaths() throws Exception {
		sumLikelihoodSets("test/input/validPathsModelAllPathsBasic.xml", train1, test1);
	}
	
	/** Tests all 2-state combinations with a restricted model and onlyValidTransitions = false */
	public void testSumPathsAllPathsCached() throws Exception {
		sumLikelihoodSets("test/input/validPathsModelAllPaths.xml", train1, test1);
	}

	/** Tests that the sum of the likelihoods of the valid paths is 1 */
	public void testSumPaths() throws Exception {
		List<? extends TrainingSequence<?>> train = StringInput.prepareData("010101\nATAGAG\n101010\nAGAGAG\n");
		String[] tests = new String[] {
				"101010101010101010101010\nATGACAGTAGACATGACAGTAGAC\n010101010101010101010101\nAATGACTGACACAATGACTGACAC\n",
				"101010101010\nATGACAGTAGAC\n010101010101\nAATGACTGACAC\n",
				"10\nAT\n01\nAA\n",
				"1\nA\n0\nT\n"
		};
		sumLikelihoodSets("test/input/validPathsModel.xml", train, tests);
	}

	/** Tests that the sum of the likelihoods of the valid paths is 1, uses a model that allows 0-0 */
	public void testSumPathsMoreComplicated() throws Exception {
		List<? extends TrainingSequence<?>> train = StringInput.prepareData("000001\nATAGAG\n101010\nAGAGAG\n");
		
		String[] tests = new String[] {
				"1\nA\n0\nT\n",
				"10\nAT\n01\nAA\n00\nTA\n",
				"000\nATG\n001\nAAC\n010\nTAG\n100\nGAG\n101\nTGG\n"
		};
		sumLikelihoodSets("test/input/validPathsModelMoreComplicated.xml", train, tests);
	}

	/* If we have an invalid path in training data we need to reject */
	public void testRejectInvalidTrainingPath() throws Exception {
		boolean failure = false;
		List<? extends TrainingSequence<?>> train = StringInput.prepareData("010101\nATAGAG\n011101\nAGAGAG\n");
		Conrad r = new Conrad("test/input/validPathsModel.xml");
		try {
			r.train(train);
		} catch (Exception e) {
			log.info(e);
			failure = true;
		}
		assertTrue(failure);
	}

	void sumLikelihoodSets(String runner, List<? extends TrainingSequence<?>> train, String[] tests) throws Exception {
		Conrad r = new Conrad(runner);
		// Test weights of 1
		r.trainFeatures(train);
		double[] weights = new double[r.getModel().getNumFeatures()];
		Arrays.fill(weights, 1.0);
		sumLikelihood(r, weights, tests);			
		
		// Test random weights
		for(int i =0; i<weights.length; ++i) {
			weights[i] = Math.random() * 10;
		}
		sumLikelihood(r, weights, tests);			
	}
	
	void sumLikelihood(Conrad r, double[] weights, String[] tests) throws Exception {
		log.debug("Sumlikelihoods test");
		for(String testString: tests) {
			// Sum a bunch of likelihoods and make sure they add up to 1.
			List<? extends TrainingSequence<?>> test = StringInput.prepareData(testString);
			double likelihood = 0.0;
			for(TrainingSequence<?> t : test) {
				CRFObjectiveFunctionGradient obj = ((StandardOptimizer)r.getOptimizer()).getObjectiveFunction();
				obj.setTrainingData(r.getModel(), Collections.singletonList(t));
				double[] dummy = new double[weights.length];
				double ll = obj.apply(weights, dummy);
				double l = Math.exp(ll*t.length());
				log.debug("Likelihood: "+l+" for "+t);
				likelihood += l;
			}
			assertEquals(1.0,likelihood,0.000001);
		}
	}
	

}
