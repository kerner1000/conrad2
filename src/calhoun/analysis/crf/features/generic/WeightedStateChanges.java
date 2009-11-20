package calhoun.analysis.crf.features.generic;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.AbstractFeatureManager;
import calhoun.analysis.crf.CacheStrategySpec;
import calhoun.analysis.crf.FeatureList;
import calhoun.analysis.crf.FeatureManagerEdge;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.CacheStrategySpec.CacheStrategy;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.util.Assert;

/** learns transition probabilities from the data and then creates a single feature for all edge transitions out of a given state.
 * <p>
 * <b>Notes:</b>
 * <ul>
 * <li> This differs from WeightedEdgeFeatures because transitions that stay in the same state are not included.  This is useful for explicit length models when
 * transitions that stay in the same state are modeled separately.   
 * <li> Returns the log probability of the transitions.  
 * <li> Uses a CONSTANT cache strategy since the values returns are independent of position.
 * <li> When learning, initializes each edge count with a pseudocount of 1.
 * </ul>
 */
public class WeightedStateChanges extends AbstractFeatureManager<Object> implements FeatureManagerEdge<Object> {
	private static final long serialVersionUID = 8477631359065280630L;
	private static final Log log = LogFactory.getLog(WeightedStateChanges.class);
	boolean debug = log.isDebugEnabled();
	
	int startIx;
	ModelManager manager;
	float[][] transitions;

	@Override
	public CacheStrategySpec getCacheStrategy() {
		return new CacheStrategySpec(CacheStrategy.CONSTANT);
	}

	public String getFeatureName(int featureIndex) {
		Assert.a(featureIndex == startIx, "Invalid feature index: ", featureIndex, ". Must be ", startIx);
		return "WeightedEdges";
	}

	public int getNumFeatures() {
		// The is a single feature that is the log transition probabilities from the data. 
		return 1;
	}

	public void evaluateEdge(InputSequence<?> seq, int pos, int prevState, int state, FeatureList result) {
		if (prevState != state) {
			result.addFeature(startIx, transitions[prevState][state]);
		}
	}

	public void train(int startingIndex, ModelManager modelInfo, List<? extends TrainingSequence<?>> data) {
		startIx = startingIndex;
		manager = modelInfo;
		int nStates = manager.getNumStates();
		
		// Count transitions from the training data
		transitions = new float[nStates][nStates];
		for (int j=0; j<nStates; j++) {
			for (int k=0; k<nStates; k++) {
				if (j!=k) {
				transitions[j][k] = (float) 1.0; // pseudocounts
				} else {
					transitions[j][k] = (float) 0.0;
				}
			}
		}
		
		//DoubleMatrix2D transitions  = new DenseDoubleMatrix2D(nStates, nStates);
		for(TrainingSequence<?> seq : data) {
			// Start at 1 because there is no transition for the first element of the sequence.
			for(int pos = 1; pos < seq.length(); ++pos) {
				int start = seq.getY(pos-1);
				int end = seq.getY(pos);
				if (start != end) {
					transitions[start][end] += (float) 1.0;
				}
			}
		}

		for (int j=0; j<nStates; j++) {
			float rowtotal = (float) 0.0;
			for (int k=0; k<nStates; k++) {
				rowtotal += transitions[j][k];
			}
			for (int k=0; k<nStates; k++) {
				if (j!= k) {
					transitions[j][k] = (float) Math.log(transitions[j][k] / rowtotal);
				}
			}
		}
	}
}
