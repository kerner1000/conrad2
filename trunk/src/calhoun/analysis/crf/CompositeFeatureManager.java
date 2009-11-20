package calhoun.analysis.crf;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import calhoun.analysis.crf.CacheStrategySpec.CacheStrategy;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.util.Assert;

/** a feature manager that combines feature types together. */
public class CompositeFeatureManager extends AbstractFeatureManager implements FeatureManagerEdge, FeatureManagerNode, FeatureManagerEdgeExplicitLength, FeatureManagerNodeExplicitLength {
	private static final long serialVersionUID = 5061912595256694049L;

	protected List<FeatureManager> allFeatureTypes = new ArrayList<FeatureManager>();

	protected List<FeatureManagerNode> nodeFeatureTypes = new ArrayList<FeatureManagerNode>();
	protected List<FeatureManagerEdge> edgeFeatureTypes = new ArrayList<FeatureManagerEdge>();
	protected List<FeatureManagerNodeExplicitLength> explicitLengthNodeFeatureTypes = new ArrayList<FeatureManagerNodeExplicitLength>();
	protected List<FeatureManagerEdgeExplicitLength> explicitLengthEdgeFeatureTypes = new ArrayList<FeatureManagerEdgeExplicitLength>();
	int[] startIndexes = null;
	
	int startIx;
	int totalFeatures;

	public 	List<FeatureManager> getComponentFeatures() {
		return allFeatureTypes;
	}
	
	public 	void setComponentFeatures(List<FeatureManager> components) {
		for ( FeatureManager component : components ) {
			addFeatureManager(component);
		}
	}
	
	public void addFeatureManager(FeatureManager fm) {
		addFeatureManager(null, null, fm);
	}
	
	public void addFeatureManager(String name, String inputParams, FeatureManager fm) {
		Assert.a(startIndexes == null, "Attempted to add a new FeatureManager after training.");

		if(name != null)
			fm.setInputComponent(name);
		allFeatureTypes.add(fm);

		// Add each feature type into the right list for evalution
		if(fm instanceof FeatureManagerNode) {
			nodeFeatureTypes.add((FeatureManagerNode) fm);
		}
		if(fm instanceof FeatureManagerEdge) {
			edgeFeatureTypes.add((FeatureManagerEdge) fm);
		}
		if(fm instanceof FeatureManagerNodeExplicitLength) {
			explicitLengthNodeFeatureTypes.add((FeatureManagerNodeExplicitLength) fm);
		}
		if(fm instanceof FeatureManagerEdgeExplicitLength) {
			explicitLengthEdgeFeatureTypes.add((FeatureManagerEdgeExplicitLength) fm);
		}
	}

	public int getNumFeatures() {
		Assert.a(startIndexes != null, "Attempted to get number of features before training.");
		return totalFeatures;
	}
	
	public String getFeatureName(int featureIndex) {
		int index = Arrays.binarySearch(startIndexes, featureIndex);
		if(index < 0) {
			index = -index-2;
		}
		// Feature managers may have no features (just constraints).
		while(allFeatureTypes.get(index).getNumFeatures() == 0) {
			index += 1;
		}
		return allFeatureTypes.get(index).getFeatureName(featureIndex);
	}

	public int getFeatureOffset(int featureIndex)
	{
		int index = Arrays.binarySearch(startIndexes, featureIndex);
		if (index < 0) {
			index = -index-2;
		}
		// Feature managers may have no features (just constraints).
		while(allFeatureTypes.get(index).getNumFeatures() == 0) {
			index += 1;
		}
		return (startIndexes[index]);
	}
	
	protected static class ComponentList extends ArrayList<InputSequence<?>> {
		private static final long serialVersionUID = 2954775229926328434L;
		String name;
		public ComponentList(List<TrainingSequence<?>> start, String name) {
			super(start);
			this.name = name;
		}
		@Override
		public InputSequence<?> get(int i) {
			return transform((TrainingSequence<?>) super.get(i));
		}
		public InputSequence<?> transform(TrainingSequence<?> in) {
			return in.getComponent(name);
		}
	}

	public void train(int startingIndex, ModelManager modelInfo, List data) {
		Assert.a(allFeatureTypes.size() > 0, "No features types have been assigned.");
		Assert.a(startIndexes == null, "FeatureManager has already been trained.");
		startIx = startingIndex;

		// Train each of the individual FeatureManagers and calculate offsets
		startIndexes = new int[allFeatureTypes.size()];
		totalFeatures = 0;
		for(int i = 0; i<startIndexes.length; ++i) {
			startIndexes[i] = totalFeatures + startIx;
			FeatureManager fm = allFeatureTypes.get(i);
			List compData = fm.getInputComponent() == null ? data : new ComponentList(data, fm.getInputComponent());
			fm.train(totalFeatures, modelInfo, compData);
			totalFeatures += fm.getNumFeatures();
		}
	}

	public void evaluateNode(InputSequence seq, int pos, int state, FeatureList result) {
		for(FeatureManagerNode fm : nodeFeatureTypes) {
			InputSequence componentSeq = fm.getInputComponent() == null ? seq : seq.getComponent(fm.getInputComponent());
			fm.evaluateNode(componentSeq, pos, state, result);
			if(!result.isValid())
				break;
		}
	}

	public void evaluateEdge(InputSequence seq, int pos, int prevState, int state, FeatureList result) {
		for(FeatureManagerEdge fm : edgeFeatureTypes) {
			InputSequence componentSeq = fm.getInputComponent() == null ? seq : seq.getComponent(fm.getInputComponent());
			fm.evaluateEdge(componentSeq, pos, prevState, state, result);
			if(!result.isValid()) {
				break;
			}
		}
	}

	public void evaluateNodeLength(InputSequence seq, int pos, int length, int state, FeatureList result) {
		Assert.a(length>0);
		for(FeatureManagerNodeExplicitLength fm : explicitLengthNodeFeatureTypes) {
			InputSequence componentSeq = fm.getInputComponent() == null ? seq : seq.getComponent(fm.getInputComponent());
			fm.evaluateNodeLength(componentSeq, pos, length, state, result);
			if(!result.isValid())
				break;
		}
	}

	public void evaluateEdgeLength(InputSequence seq, int pos, int length, int prevState, int state, FeatureList result) {
		Assert.a(length>0);
		for(FeatureManagerEdgeExplicitLength fm : explicitLengthEdgeFeatureTypes) {
			InputSequence componentSeq = fm.getInputComponent() == null ? seq : seq.getComponent(fm.getInputComponent());
			fm.evaluateEdgeLength(componentSeq, pos, length, prevState, state, result);
			if(!result.isValid())
				break;
		}
	}

	@Override
	public CacheStrategySpec getCacheStrategy() {
		return new CacheStrategySpec(CacheStrategy.COMPOSITE);
	}
}
