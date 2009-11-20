package calhoun.analysis.crf.test;

import java.util.List;

import calhoun.analysis.crf.AbstractFeatureManager;
import calhoun.analysis.crf.CacheStrategySpec;
import calhoun.analysis.crf.FeatureList;
import calhoun.analysis.crf.FeatureManager;
import calhoun.analysis.crf.FeatureManagerEdge;
import calhoun.analysis.crf.FeatureManagerNode;
import calhoun.analysis.crf.FeatureManagerNodeExplicitLength;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.CacheStrategySpec.CacheStrategy;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.analysis.crf.io.TrainingSequence;

public class TestFeatures {
	public static abstract class TestFeature extends AbstractFeatureManager<Character> implements FeatureManager<Character> {

		int startIx;
		
		public String getFeatureName(int featureIndex) {
			return getClass().getName();
		}

		public int getNumFeatures() {
			return 1;
		}

		public void train(int startingIndex, ModelManager modelInfo, List<? extends TrainingSequence<? extends Character>> data) {
			startIx = startingIndex;
		}
	}

	/* Returns log(.8) for G/C in state 1 or A/T in state 0, log(.2) otherwise */
	public static class EmissionFeature extends TestFeature implements FeatureManagerNode<Character> {
		private static final long serialVersionUID = 8039479168741225007L;

		public void evaluateNode(InputSequence<? extends Character> seq, int pos, int state, FeatureList result) {
			Character c = seq.getX(pos);
			boolean match = state == 0;
			if(c.charValue() == 'G' || c.charValue() == 'C') {
				match = state == 1;
			}
			result.addFeature(startIx, match ? Math.log(.8) : Math.log(.2));
		}
		@Override
		public CacheStrategySpec getCacheStrategy() {
			return new CacheStrategySpec(CacheStrategy.UNSPECIFIED);
		}
	}

	/* Returns 0 for state 0 transitions.  Returns guassian with mena of 6 and std dev of 1. for state 1 */
	public static class GaussianLengthFeature extends TestFeature implements FeatureManagerNodeExplicitLength<Character> {
		private static final long serialVersionUID = 6050417482057409153L;

		public void evaluateNodeLength(InputSequence<? extends Character> seq, int pos, int length, int state, FeatureList result) {
			double val = 0.0;
			if(state == 1) {
				val = Math.pow(length - 6, 2);
				result.addFeature(startIx, val);
			}
		}
		@Override
		public CacheStrategySpec getCacheStrategy() {
			return new CacheStrategySpec(CacheStrategy.UNSPECIFIED);
		}
	}

	/* Has a value of 1/2 the edge feature value for staying the same transition. */
	public static class ExplicitHalfExponentialLengthFeature extends TestFeature implements FeatureManagerNodeExplicitLength<Character>, FeatureManagerEdge<Character> {
		private static final long serialVersionUID = 6050417482057409153L;

		public void evaluateNodeLength(InputSequence<? extends Character> seq, int pos, int length, int state, FeatureList result) {
			if(state == 1) {
				float val = (float)((length-1)*Math.log(.8)/2.0f);
				result.addFeature(startIx, val);
			}
		}

		float[] vals = new float[] {.8f, .2f, .2f, (float) Math.sqrt(.8)};
		
		public void evaluateEdge(InputSequence<? extends Character> seq, int pos, int prevState, int state, FeatureList result) {
			result.addFeature(startIx, Math.log(vals[prevState + 2*state]));
		}
		@Override
		public CacheStrategySpec getCacheStrategy() {
			return new CacheStrategySpec(CacheStrategy.UNSPECIFIED);
		}
	}

	public static class ExponentialLengthFeature extends TestFeature implements FeatureManagerEdge<Character> {
		private static final long serialVersionUID = 290546936098017942L;

		float[] vals = new float[] {.8f, .2f, .2f, .8f};
		
		public void evaluateEdge(InputSequence<? extends Character> seq, int pos, int prevState, int state, FeatureList result) {
			result.addFeature(startIx, Math.log(vals[prevState + 2*state]));
		}
		@Override
		public CacheStrategySpec getCacheStrategy() {
			return new CacheStrategySpec(CacheStrategy.UNSPECIFIED);
		}
	}

	public static class HalfExponentialLengthFeature extends TestFeature implements FeatureManagerEdge<Character> {
		private static final long serialVersionUID = 290546936098017942L;

		float[] vals = new float[] {.8f, .2f, .2f, (float) Math.sqrt(.8)};
		
		public void evaluateEdge(InputSequence<? extends Character> seq, int pos, int prevState, int state, FeatureList result) {
			result.addFeature(startIx, Math.log(vals[prevState + 2*state]));
		}
		@Override
		public CacheStrategySpec getCacheStrategy() {
			return new CacheStrategySpec(CacheStrategy.UNSPECIFIED);
		}
	}
}
