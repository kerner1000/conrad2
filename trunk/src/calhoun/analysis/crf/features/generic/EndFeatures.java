package calhoun.analysis.crf.features.generic;

import java.util.Arrays;
import java.util.List;

import calhoun.analysis.crf.AbstractFeatureManager;
import calhoun.analysis.crf.FeatureList;
import calhoun.analysis.crf.FeatureManagerNode;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.util.Assert;

/** indicator functions that evaluate to true for a selected set of start states at the first position in the sequence. */
public class EndFeatures extends AbstractFeatureManager<Object> implements FeatureManagerNode<Object> {
	private static final long serialVersionUID = 5206582215565363768L;
	int startIx;
	int[] endStates;
	String[] names;
	
	int[] configStates = new int[] {0, 1};

	public void setStates(int[] config) {
		configStates = config;
	}
	
	public String getFeatureName(int featureIndex) {
		Assert.a(featureIndex - startIx < names.length, "Invalid feature index");
		return names[featureIndex - startIx];
	}

	public int getNumFeatures() {
		return names.length;
	}

	public void evaluateNode(InputSequence<?> seq, int pos, int state, FeatureList result) {
		if(pos == seq.length()-1) {
			int index = endStates[state];
			if(index != -1) {
				result.addFeature(index, 1);
			}
		}
	}

	/** Edge features don't train based on the data.  Just set up based on the model. */
	public void train(int startingIndex, ModelManager modelInfo, List data) {
		startIx = startingIndex;
		endStates = new int[modelInfo.getNumStates()];
		Arrays.fill(endStates, -1);
		names = new String[configStates.length];
		for(int i=0; i<configStates.length; ++i) {
			endStates[configStates[i]] = startIx+i;
			names[i] = "End."+modelInfo.getStateName(configStates[i]);
		}
	}
}
