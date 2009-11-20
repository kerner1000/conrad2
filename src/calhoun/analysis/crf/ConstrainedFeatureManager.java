package calhoun.analysis.crf;

import java.util.List;

import org.apache.commons.lang.StringUtils;

import calhoun.util.Assert;

/** a feature manager that combines it's composite feature types together into a single feature. */
public class ConstrainedFeatureManager extends CompositeFeatureManager {
	private static final long serialVersionUID = 5061912595256694050L;

	/** Composite feature managers only return 1 feature. */
	@Override
	public int getNumFeatures() {
		return 1;
	}

	/** Return a concatenated name */
	@Override
	public String getFeatureName(int featureIndex) {
		String[] names = new String[allFeatureTypes.size()];
		for(int i=0; i<names.length; ++i) {
			names[i] = allFeatureTypes.get(i).getFeatureName(featureIndex);
		}
		return "Composite: "+StringUtils.join(names, ",");
	}

	@Override
	public void train(int startingIndex, ModelManager modelInfo, List data) {
		Assert.a(allFeatureTypes.size() > 0, "No features types have been assigned.");
		Assert.a(startIndexes == null, "FeatureManager has already been trained.");
		startIx = startingIndex;

		// Train each of the individual FeatureManagers and calculate offsets
		for(int i = 0; i<allFeatureTypes.size(); ++i) {
			FeatureManager fm = allFeatureTypes.get(i);
			List compData = fm.getInputComponent() == null ? data : new ComponentList(data, fm.getInputComponent());
			fm.train(startIx, modelInfo, compData);
			Assert.a(fm.getNumFeatures()==1, "Constrained FeatureManagers must all have 1 feature.  ",fm, " had ", fm.getNumFeatures());
		}
	}
}
