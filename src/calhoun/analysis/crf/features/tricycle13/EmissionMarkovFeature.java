package calhoun.analysis.crf.features.tricycle13;

import java.io.Serializable;
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
import calhoun.analysis.crf.features.supporting.MarkovPredictorLogprob;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.analysis.crf.io.TrainingSequence;

public class EmissionMarkovFeature extends AbstractFeatureManager<Character> implements FeatureManagerNode<Character> {
	private static final long serialVersionUID = -7659288739348604129L;
	private static final Log log = LogFactory.getLog(EmissionMarkovFeature.class);
	boolean debug = log.isDebugEnabled();
	
	int startIx;  // The index of the first feature managed by this FeatureManager
	ModelManager model;
	
	// Following block are things that depend explicitly on and are calculated directly
	//  from geometry, included only for convenience. 
	MarkovPredictorLogprob predictorlp;
	boolean tieFlag = false;

	public static class MarkovHistory implements Serializable {
		private static final long serialVersionUID = 8731309130784681776L;

		List<List<Node>> states;

		public List<List<Node>> getStates() {
			return states;
		}

		public void setStates(List<List<Node>> history) {
			this.states = history;
		}
		
		public List<int[]> convert() {
			List<int[]> historyArray= new ArrayList<int[]>();
			for(List<Node> nodeList : states) {
				int[] historyEntry = new int[nodeList.size()];
				historyArray.add(historyEntry);
				for(int i=0; i<historyEntry.length; ++i) {
					historyEntry[i] = nodeList.get(i).getIndex();
				}
			}
			return historyArray;
		}
	}
	
	public EmissionMarkovFeature() {
	}

	public void setHistory(MarkovHistory markovHistory) {
		this.predictorlp = new MarkovPredictorLogprob(markovHistory.convert());
	}
	
	public EmissionMarkovFeature(List<int[]> history) {
		this.predictorlp = new MarkovPredictorLogprob(history);
	}

	public EmissionMarkovFeature(List<int[]> history, List<int[]> flags) {
		this.predictorlp = new MarkovPredictorLogprob(history);
		tieFlag = true;
	}

	public int getNumFeatures() {
		return tieFlag ? 1 : model.getNumStates();
	}	
	
	public String getFeatureName(int featureIndex) {
		if (tieFlag) { 
			return "TiedEmissionMarkovFeature";
		} else {
			int raw = featureIndex - startIx;		
			return "EmissionMarkov.span" + model.getStateName(raw);
		}
	}
	
	public void evaluateNode(InputSequence<? extends Character> seq, int pos, int state, FeatureList result) {
		int index = startIx + (tieFlag ? 0 : state);
		result.addFeature(index, predictorlp.logprob(state,seq,pos));
	}
	
	public void train(int startingIndex, ModelManager modelInfo, List<? extends TrainingSequence<? extends Character>> data) {
		startIx = startingIndex;
		model = modelInfo;
		
		predictorlp.train(data);
	}

	@Override
	public CacheStrategySpec getCacheStrategy() {
		return new CacheStrategySpec(CacheStrategy.DENSE);
	}
}


