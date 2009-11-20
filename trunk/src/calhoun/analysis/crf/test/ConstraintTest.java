package calhoun.analysis.crf.test;

import java.util.List;

import calhoun.analysis.crf.AbstractFeatureManager;
import calhoun.analysis.crf.CacheStrategySpec;
import calhoun.analysis.crf.Conrad;
import calhoun.analysis.crf.FeatureList;
import calhoun.analysis.crf.FeatureManagerEdge;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.CacheStrategySpec.CacheStrategy;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.util.AbstractTestCase;

public class ConstraintTest extends AbstractTestCase {
	public void testViterbiConstraints() throws Exception {
		Conrad r = new Conrad("test/input/constraintModel.xml");
		r.trainFeatures("test/input/shortTrain.tricycle13.txt");
		r.setWeights(new double[] {1.0});
		
		// Test that the feature name prints out correctly
		assertEquals("Edge Feature", r.getModel().getFeatureName(0));
		
		// Test that the constraints are observed by checking against a manually validated set.
		r.test("test/input/shortTrain.tricycle13.txt", "test/working/constraintTest.txt");
		assertFilesMatch("test/output/constraintTest.txt", "test/working/constraintTest.txt");
	}

	
	public void testConstraintsDifferentCaches() throws Exception {
		Conrad r = new Conrad("test/input/constraintModel.xml");
		r.train("test/input/genesShorterTrain.txt");
		double normWeight = r.getWeights()[0];
		
		r = new Conrad("test/input/constraintModelCache.xml");
		r.train("test/input/genesShorterTrain.txt");
		assertEquals(normWeight, r.getWeights()[0], 0.0001);
	}
	
	/** Edge class used for testing that favors changing states whenever possible.  Used to generate interesting paths for constraints */
	public static class FixedEdges extends AbstractFeatureManager implements FeatureManagerEdge {
		private static final long serialVersionUID = 5995552526733022868L;
		int startIx;
		int[] clusters = new int[] {0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4};
		public void evaluateEdge(InputSequence seq, int pos, int prevState, int state, FeatureList result) {
			result.addFeature(startIx, clusters[state] == clusters[prevState] ? Math.log(.3) : Math.log(.7));
		}

		public String getFeatureName(int featureIndex) {
			return "Edge Feature";
		}

		public int getNumFeatures() {
			return 1;
		}

		public void train(int startingIndex, ModelManager modelInfo, List data) {
			startIx = startingIndex;
		}
	
		@Override
		public CacheStrategySpec getCacheStrategy() {
			return new CacheStrategySpec(CacheStrategy.UNSPECIFIED);
		}
			
	}
	
	
}
