package calhoun.analysis.crf.features.tricycle13;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.AbstractFeatureManager;
import calhoun.analysis.crf.CacheStrategySpec;
import calhoun.analysis.crf.FeatureList;
import calhoun.analysis.crf.FeatureManagerNode;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.BeanModel.Node;
import calhoun.analysis.crf.CacheStrategySpec.CacheStrategy;
import calhoun.analysis.crf.features.supporting.phylogenetic.ColumnConditionalLogProbability;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.analysis.crf.io.MultipleAlignmentInputSequence.MultipleAlignmentColumn;

public class FelsensteinFeatures extends AbstractFeatureManager<MultipleAlignmentColumn> implements FeatureManagerNode<MultipleAlignmentColumn> {
	private static final long serialVersionUID = -7659288739348604129L;
	private static final Log log = LogFactory.getLog(FelsensteinFeatures.class);
	boolean debug = log.isDebugEnabled();
	
	/* We implement the following features:
	   feature_intronicp(y_i-1,y_i,x,i)
	     = log(pr(multiplealignmentcolumn | y_i = {intronicp}, reference seq nuceotide, evolutionary model))
	         * delta(y_i=intronicp),
	   and similar features for collections of states other than intronicp = {intron1,intron2,intron3},
	   as specified by configuration properties.
	   
	   Given the evolutionary model, including a phylogenetic tree, this probability can be
	   evaluated efficiently, even with missing data, using Felsenstein's algorithm.
	   
	   The evolutionary model must be trained using maximum likelihood.
	   For this, we will use a Nelder-Mead solver that requires function evaluations but
	   does not require gradient evaluations.
	
	   We will take the topology and the relative branch lengths of the phylogenetic tree as given
	   (perhaps in the model configuration file).  It is the overall scaling and parameters such as
	   (for an HKY model) the ratio of transitions to transversions that must be determined by maximum
	   likelihood.
	   
	   Choice looming ahead: where should the phylogenetic tree with relative branchlengths be represented?
	     a) within the Multiple alignment, ie input from file from the data
	     b) within the model configuration file
	     c) a third place altogether different
	   I think maybe choice a is best...but let's be open minded until need to choose.
	*/

	int startIx;  // The index of the first feature managed by this FeatureManager
	ModelManager model;	
	ColumnConditionalLogProbability mo;
	boolean tieFlag = false;
	List<int[]> clusterIndices;
	private int eModelNum = 0; // See ColumnConditionalLogProbability for interpetation; 0 is the default.
	
	public FelsensteinFeatures(List<int[]> clusters) throws ClassNotFoundException {
		this.clusterIndices = clusters;
	}

	public void setClusters(List<List<Node>> clusters) {
		clusterIndices = new ArrayList();
		for(List<Node> nodeList : clusters) {
			int[] cluster = new int[nodeList.size()]; 
			clusterIndices.add(cluster);
			for(int i=0; i<nodeList.size(); ++i) {
				cluster[i] = nodeList.get(i).getIndex();
			}
		}
	}
	
	public FelsensteinFeatures() { }	
	
	public FelsensteinFeatures(List<int[]> clusters, List<int[]> eModelNum) {
		this.clusterIndices = clusters;
		this.eModelNum = eModelNum.get(0)[0];
	}
	
	public FelsensteinFeatures(List<int[]> clusters, List<int[]> eModelNum, List<int[]> flags) {
		tieFlag = true;
		this.clusterIndices = clusters;
		this.eModelNum = eModelNum.get(0)[0];
	}
	
	public int getNumFeatures() {
		if (tieFlag) {
			return 1;
		} else {
			return mo.numClusters();			
		}
	}	
	
	public String getFeatureName(int featureIndex) {
		return "FelsensteinFeatures";
	}


	public void evaluateNode(InputSequence<? extends MultipleAlignmentColumn> seq, int pos, int state, FeatureList result) {
		int cl;
		if (tieFlag) {
			cl = 0;
		} else {
			cl = mo.state2cluster(state);
		}
		result.addFeature(startIx + cl, mo.condLogProb(seq,pos,state));
	}

	
	public void train(int startingIndex, ModelManager modelInfo, final List<? extends TrainingSequence<? extends MultipleAlignmentColumn>> data) {
		startIx = startingIndex;
		model = modelInfo;
	
		mo = new ColumnConditionalLogProbability(clusterIndices,eModelNum);
		
		mo.train(modelInfo,data);
	}
	@Override
	public CacheStrategySpec getCacheStrategy() {
		return new CacheStrategySpec(CacheStrategy.UNSPECIFIED);
	}
}
