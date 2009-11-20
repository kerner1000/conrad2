/**
 * 
 */
package calhoun.analysis.crf.test;

import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.AbstractFeatureManager;
import calhoun.analysis.crf.CacheStrategySpec;
import calhoun.analysis.crf.CompositeFeatureManager;
import calhoun.analysis.crf.FeatureList;
import calhoun.analysis.crf.FeatureManagerNode;
import calhoun.analysis.crf.FeatureManagerNodeBoundaries;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.CacheStrategySpec.CacheStrategy;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.util.Assert;
import calhoun.util.DenseBooleanMatrix2D;

/** Test manager that avoids assigning length based features and non-length based features the same index. 
 * Also allows testing of Dense Node features*/
class TestFeatureManager3 extends CompositeFeatureManager implements ModelManager {
	private static final long serialVersionUID = 5162614677151226890L;
	@SuppressWarnings("unused")
	private static final Log log = LogFactory.getLog(TestFeatureManager3.class);

	boolean len;
	
	public TestFeatureManager3(int type, boolean len) {
		this.len = len;
		// type 0 - 1 edge feature
		// type 1 - 2 node features
		// type 2 - 2 node features + 1 edge feature
		// type 3 - 1 node feature (length or no length)
		// If len is true, second node feature is length based.
		switch(type) {
		case 1:
			addFeatureManager(new TestNodeFeature());
		case 3:
			break;
		default:
			Assert.a(false, "Type "+type+" not 1 or 3.");
		}
		if(len)
			addFeatureManager(new DenseTestNodeLengthFeature());
		else
			addFeatureManager(new TestNodeFeature());
	}

	abstract class TestFeature extends AbstractFeatureManager implements FeatureManagerNode {
		public void evaluateNode(InputSequence seq, int pos, int state, FeatureList result) {
			double val = (pos%5+1)/6.0f;
			if(state == 1) {
				val = 1 - val;
			}
			val = Math.log(val) ;//* 0.5;
			result.addFeature(startIx, val);
			//log.info(String.format("F: %d P: %d S: %d Val: %f", startIx, pos, state, val));
		}

		int startIx;
		
		public String getFeatureName(int featureIndex) {
			Assert.a(featureIndex == startIx);
			return "DenseNode";
		}

		public void train(int startingIndex, ModelManager modelInfo, List data) {
			startIx = startingIndex;
		}

		public int getNumFeatures() {
			return 1;
		}
	}

	class TestNodeFeature extends TestFeature implements FeatureManagerNode {
		private static final long serialVersionUID = 5966577419490935155L;

		@Override
		public CacheStrategySpec getCacheStrategy() {
			CacheStrategySpec css = new CacheStrategySpec(CacheStrategy.DENSE);
			return css;
		}
	}

	class DenseTestNodeLengthFeature extends TestFeature implements FeatureManagerNodeBoundaries {
		private static final long serialVersionUID = -8448250464759172704L;

		@Override
		public CacheStrategySpec getCacheStrategy() {
			CacheStrategySpec css = new CacheStrategySpec(CacheStrategy.DENSE_NODE_BOUNDARY);
			CacheStrategySpec.DenseBoundaryCachingDetails details = new CacheStrategySpec.DenseBoundaryCachingDetails(2); // we will use 9 tables
			
			// state, table, featureIx, rightPad, leftPad
			// If you want to debug and verify equality with the Markov case, then set all the pads to zero as below:
			
			details.add(0,0,startIx,0,0);
			details.add(1,1,startIx,0,0);

			details.check();
					
			css.details = details;
			
			return css;
		}
	}


	public int getNumStates() {
		return 2;
	}

	public String getStateName(int state) {
		return "testState"+state;
	}

	public int getStateIndex(String name) {
		Assert.a(false);
		return 0;
	}

	public DenseBooleanMatrix2D getLegalTransitions() {
		if(len) {
			DenseBooleanMatrix2D trans = new DenseBooleanMatrix2D(2, 2);
			trans.setQuick(0, 1, true);
			trans.setQuick(1, 0, true);
			return trans;
		}
		return null;
	}

}