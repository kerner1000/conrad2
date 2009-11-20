package calhoun.analysis.crf.features.generic;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.AbstractFeatureManager;
import calhoun.analysis.crf.FeatureList;
import calhoun.analysis.crf.FeatureManagerEdge;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.util.Assert;
import calhoun.util.DenseBooleanMatrix2D;
import calhoun.util.DenseIntMatrix2D;

/** a set of indicator features for the valid transitions in the model.
 * <p>
 * <b>Notes:</b>
 * <ul>
 * <li> Each indicator returns 1 when evaluated on it's transition, otherwise 0.
 * <li> Exactly one feature will return a non-zero value for every valid edge transition. 
 * </ul>
 */
public class IndicatorEdges extends AbstractFeatureManager implements FeatureManagerEdge {
	private static final long serialVersionUID = -2633500053477439285L;
	private static final Log log = LogFactory.getLog(IndicatorEdges.class);
	boolean debug = log.isDebugEnabled();
	
	int startIx;
	DenseIntMatrix2D transitions;
	List<String> names;
	
	public String getFeatureName(int featureIndex) {
		Assert.a(featureIndex - startIx < names.size(), "Invalid feature index");
		return names.get(featureIndex - startIx);
	}

	public int getNumFeatures() {
		return names.size();
	}

	public void evaluateEdge(InputSequence seq, int pos, int prevState, int state, FeatureList result) {
		int index = transitions.getQuick(prevState, state);
		if(index != -1) {
			result.addFeature(index, 1);
		}
	}

	/** Edge features don't train based on the data.  Just set up based on the model. */
	public void train(int startingIndex, ModelManager modelInfo, List data) {
		startIx = startingIndex;
		int nStates = modelInfo.getNumStates();
		transitions = new DenseIntMatrix2D(nStates, nStates);
		transitions.assign(-1);
		names = new ArrayList<String>();
		
		int n = 0;
		DenseBooleanMatrix2D trans = modelInfo.getLegalTransitions();
		for(int i=0; i<nStates; ++i) {
			for(int j=0; j<nStates; ++j) {
				if(trans.getQuick(i, j)) {
					transitions.setQuick(i, j, startIx + n);
					names.add("Edge."+modelInfo.getStateName(i)+"-"+modelInfo.getStateName(j));
					++n;
				}
			}
		}
	}
}
