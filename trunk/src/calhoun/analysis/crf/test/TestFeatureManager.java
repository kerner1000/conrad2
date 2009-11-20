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

class TestFeatureManager extends AbstractFeatureManager implements ModelManager {
	private static final long serialVersionUID = 5316261467715122688L;
	@SuppressWarnings("unused")
	private static final Log log = LogFactory.getLog(TestFeatureManager.class);

	int type;
	boolean lengthModel;
	boolean lengthFeatures;
	boolean edgeLengthFeatures = false;	// Not supported right now.
	
	public TestFeatureManager(int type, int length) {
		this.type=type;
		this.lengthModel=length > 0;
		this.lengthFeatures=length == 2;
	}

	public TestFeatureManager(int type) {
		this(type, 0);
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
		return type == 2 ? 2:1;
	}

	public void train(int startingIndex, ModelManager modelInfo, List data) {
	}

	public DenseBooleanMatrix2D getLegalTransitions() {
		if(lengthModel) {
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
		if(type == 1 || type >= 2){
			// Edge features for 2/3 stay in state, 1/3 leave.
			result.addFeature(type == 2 ? 1:0, (float) Math.log(prevState == state ? 2/3.0 : 1/3.0) * (edgeLengthFeatures?0.5:1));
		}
	}

	public void evaluateNode(InputSequence seq, int pos, int state, FeatureList result) {
		if(type == 0 || type >= 2) {
			result.addFeature(0, (float) Math.log(state==0 ? (1.0f/3.0f) : (2.0f/3.0f)) * (lengthFeatures?0.5:1));
		}
	}

	/* In the length model, the regular node feature counts half and the length node feature counts half.
	 * The regular edge feature counts half for both the state change transitions and the self-transitions.
	 * The length edge feature count half of the state change plus half of the same state transitions.
	 * The net result is a model with the same behavior as the non-length features.
	 */
	public void evaluateEdgeLength(InputSequence seq, int pos, int featLength, int prevState, int state, FeatureList result) {
		if(edgeLengthFeatures && (type == 1 || type >= 2)) {
			// Edge features for 2/3 stay in state, 1/3 leave.
			double val = .5*Math.log(prevState == state ? 2/3.0 : 1/3.0) + .5*(featLength-1)*Math.log(2/3.0);
			result.addFeature(type == 2 ? 1:0, (float) val);
			//log.info(String.format("Edge: %d %d %d %d %f", pos, featLength, prevState, state, val));
		}
	}

	public void evaluateNodeLength(InputSequence seq, int pos, int featLength, int state, FeatureList result) {
		if(edgeLengthFeatures && (type == 1 || type >= 2) && featLength == pos + 1) {
			// Handle edge lengths for the features that start at the beginning
			result.addFeature(type == 2 ? 1:0, (float) .5*(featLength-1)*Math.log(2/3.0));
			//log.info(String.format("Node: %d %d %d %d", pos, featLength, -1, state));
		}
		if(lengthFeatures && (type == 0 || type >= 2)) {
			double val = .5*Math.log(state==0 ? 1.0f/3.0f : 2.0f/3.0f) * featLength;
			result.addFeature(0, (float) val);
			//log.info(String.format("Node: %d %d %d %d %f", pos, featLength, -1, state, val));
		}
	}

	public void buildFeatureCache(String directory, List sequence) {		
	}
	@Override
	public CacheStrategySpec getCacheStrategy() {
		return new CacheStrategySpec(CacheStrategy.UNSPECIFIED);
	}
}