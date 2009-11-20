/**
 * 
 */
package calhoun.analysis.crf.test;

import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.AbstractFeatureManager;
import calhoun.analysis.crf.CacheStrategySpec;
import calhoun.analysis.crf.FeatureList;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.CacheStrategySpec.CacheStrategy;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.util.Assert;
import calhoun.util.DenseBooleanMatrix2D;

/** Test manager that avoids assigning length based features and non-length based features the same index. */
class TestFeatureManager2 extends AbstractFeatureManager implements ModelManager {
	private static final long serialVersionUID = 5316261467715122688L;
	@SuppressWarnings("unused")
	private static final Log log = LogFactory.getLog(TestFeatureManager2.class);

	int edgeIx = -1;
	int node1Ix = -1;
	int node2Ix = -1;
	int nodeLenIx = -1;
	int nFeatures;
	
	public TestFeatureManager2(int type, boolean len) {
		// type 0 - 1 edge feature
		// type 1 - 2 node features
		// type 2 - 2 node features + 1 edge feature
		// type 3 - 1 node feature (length or no length)
		// If len is true, second node feature is length based.
		switch(type) {
		case 0:
			nFeatures = 1;
			edgeIx = 0;
			break;
		case 1:
			nFeatures = 2;
			node1Ix = 0;
			if(len)
				nodeLenIx = 1;
			else
				node2Ix = 1;
			break;
		case 2:
			nFeatures = 3;
			edgeIx = 0;
			node1Ix = 1;
			if(len)
				nodeLenIx = 2;
			else
				node2Ix = 2;
			break;
		case 3:
			nFeatures = 1;
			if(len)
				nodeLenIx = 0;
			else
				node2Ix = 0;
			break;
		default:
			Assert.a(false, "Type "+type+" not 0-2.");
		}
	}

	public int getNumStates() {
		return 2;
	}

	public String getStateName(int state) {
		return "testState"+state;
	}

	public String getFeatureName(int featureIndex) {
		return "testFeature"+featureIndex;
	}

	public int getNumFeatures() {
		return nFeatures;
	}

	public void train(int startingIndex, ModelManager modelInfo, List data) {
	}

	public DenseBooleanMatrix2D getLegalTransitions() {
		if(nodeLenIx != -1) {
			DenseBooleanMatrix2D trans = new DenseBooleanMatrix2D(2, 2);
			trans.setQuick(0, 1, true);
			trans.setQuick(1, 0, true);
			return trans;
		}
		return null;
	}

	public int getStateIndex(String name) {
		Assert.a(false);
		return 0;
	}

	public void evaluateEdge(InputSequence seq, int pos, int prevState, int state, FeatureList result) {
		if(edgeIx != -1){
			// Edge features for 2/3 stay in state, 1/3 leave.
			result.addFeature(edgeIx, (float) Math.log(prevState == state ? 2/3.0 : 1/3.0));
		}
	}

	public void evaluateNode(InputSequence seq, int pos, int state, FeatureList result) {
		double val = Math.log(state==0 ? (1.0f/3.0f) : (2.0f/3.0f)) * 0.5;		
		if(node1Ix != -1)
			result.addFeature(node1Ix, val);
		if(node2Ix != -1)
			result.addFeature(node2Ix, val);
		//log.info(String.format("Node: %d %d %f", pos, state, val));
	}

	/* In the length model, the regular node feature counts half and the length node feature counts half.
	 * The regular edge feature counts half for both the state change transitions and the self-transitions.
	 * The length edge feature count half of the state change plus half of the same state transitions.
	 * The net result is a model with the same behavior as the non-length features.
	 */
	public void evaluateEdgeLength(InputSequence seq, int pos, int featLength, int prevState, int state, FeatureList result) {
//		if(edgeLengthFeatures && (type == 1 || type >= 2)) {
//			// Edge features for 2/3 stay in state, 1/3 leave.
//			double val = .5*Math.log(prevState == state ? 2/3.0 : 1/3.0) + .5*(featLength-1)*Math.log(2/3.0);
//			result.addFeature(type == 2 ? 1:0, (float) val);
//			//log.info(String.format("Edge: %d %d %d %d %f", pos, featLength, prevState, state, val));
//		}
	}

	public void evaluateNodeLength(InputSequence seq, int pos, int featLength, int state, FeatureList result) {
		if(nodeLenIx != -1) {
			double val = .5*Math.log(state==0 ? 1.0f/3.0f : 2.0f/3.0f) * featLength;
			result.addFeature(nodeLenIx, (float) val);
			//log.info(String.format("Node: %d %d %d %f", pos, featLength, state, val));
		}
	}

	@Override
	public CacheStrategySpec getCacheStrategy() {
		return new CacheStrategySpec(CacheStrategy.UNSPECIFIED);
	}
}