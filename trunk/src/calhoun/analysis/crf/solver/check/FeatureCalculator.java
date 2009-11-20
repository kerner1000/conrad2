/**
 * 
 */
package calhoun.analysis.crf.solver.check;

import java.util.Arrays;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.util.Assert;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;

public class FeatureCalculator {
	private static final Log log = LogFactory.getLog(FeatureCalculator.class);
	boolean debug = log.isDebugEnabled();;
	
	public ArrayFeatureList result;
	ModelManager manager;
	double[] lambda;
	double weightedFeatureSum;
	double[] featureSums;
	int nFeatures;
	int nStates;
	TransitionInfo transitions;
	
	/** If featureSums is true, summed values are tracked for each feature.  The caller must be sure to make sure
	 * call the calculator for each node and edge exactly once.
	 * @param manager
	 * @param lambda
	 */ 
	public FeatureCalculator(ModelManager manager, double[] lambda, TransitionInfo transitions) {
		this.manager = manager;
		this.lambda = lambda;

		nFeatures = manager.getNumFeatures();
		nStates = manager.getNumStates();
		result = new ArrayFeatureList(manager);
		this.transitions = transitions;
	}

	/** Returns the current expectation */
	public double getWeightedFeatureSum() {
		return weightedFeatureSum;
	}
	
	public double[] getFeatureSums() {
		return featureSums;
	}

	/** Resets the expectation to 0 and sets the compute expectation flag */
	public void resetFeatureSums() {
		weightedFeatureSum = 0.0;
		if(featureSums == null)
			featureSums = new double[nFeatures];
		else
			Arrays.fill(featureSums, 0);
	}
	
	/** Evaluate all features for a particular edge */
	public double calcEdgeValue(InputSequence seq, int pos, int prevState, int state) {
		if(debug) {
			log.debug(String.format("Edge - pos %d prevState %d state %d", pos, prevState, state));
		}
		if(transitions.transitionIndex.getQuick(prevState, state) == -1) {
			return Double.NEGATIVE_INFINITY;
		}
		result.evaluateEdge(seq, pos, prevState, state);
		boolean updateSum = checkForUpdate(seq, pos, prevState, state);
		return calcRet(updateSum);
	}
	
	public boolean isValidTransition(int prevState, int state) {
		if (transitions.transitionIndex.getQuick(prevState, state) == -1)
			return false;
		return true;
	}

	/** Evaluate all features for a particular node */
	public double calcNodeValue(InputSequence seq, int pos, int state) {
		if(debug) {
			log.debug(String.format("Node - pos %d state %d", pos, state));
		}
		result.evaluateNode(seq, pos, state);
		boolean updateSum = checkForUpdate(seq, pos, -1, state);
		return calcRet(updateSum);
	}

	/** Update a transition matrix for transitions into and nodes at a given pos */
	public void computeMi(InputSequence seq, int pos, DoubleMatrix2D mi, DoubleMatrix1D ri) {
		for(int current=0; current < nStates; ++current) {
			double nodeVal = calcNodeValue(seq, pos, current);
			if(ri != null) {
				ri.setQuick(current, nodeVal);
			}
			if(pos > 0) {
				for(int prev=0; prev < nStates; ++prev) {
					mi.setQuick(prev, current, nodeVal + calcEdgeValue(seq, pos, prev, current));
				}
			}
		}
	}
	
	/** Update a transition matrix for transitions into and nodes at a given pos */
	public void computeSparseMi(InputSequence seq, int pos, double[] mi, double[] ri) {
		double nodeVal = 0.0;
		int currentState = -1;
		for(int current :  transitions.orderedPotentials) {
			if(current < nStates) {
				nodeVal = calcNodeValue(seq, pos, current);
				currentState = current;
				if(ri != null) {
					ri[current] = nodeVal;
				}
			}
			else if(pos > 0) {
				mi[current - nStates] = nodeVal + calcEdgeValue(seq, pos, transitions.transitionFrom[current-nStates], currentState);
			}
		}
	}
	
	public double calcNodeLengthValue(InputSequence seq, int pos, int len, int state) {
		result.evaluateNodeLength(seq, pos, len, state);
		//boolean updateSum = checkForUpdate(seq, pos, -1, state);
		return calcRet(false);
	}

	public double calcEdgeLengthValue(InputSequence seq, int pos, int len, int prevState, int state) {
		result.evaluateEdgeLength(seq, pos, len, prevState, state);
		//boolean updateSum = checkForUpdate(seq, pos, -1, state);
		return calcRet(false);
	}

	/** Checks if an expectation should be computed and if the training data matches the curretn node or edge */
	boolean checkForUpdate(InputSequence seq, int pos, int previousState, int state) {
		if(featureSums != null) {
			TrainingSequence train = (TrainingSequence) seq;
			boolean prevMatches = pos == 0 || previousState == -1 || previousState == train.getY(pos-1); 
			return prevMatches && (state == train.getY(pos));
		}
		return false;
	}

	/** Compute the weighted sum of the features values */
	public double calcRet(boolean updateSum) {
		if(!result.valid) {
			return Double.NEGATIVE_INFINITY;
		}
		double ret = 0.0;
		int count = result.currentSize; //result.size();
		int[] indices = result.indices; //result.getIndices();
		double[] vals = result.values; //result.getValues();
		for(int i = 0; i<count; ++i) {
			int index = indices[i];
			double val = vals[i] * lambda[index];
			ret += val;
			if(updateSum) {
				featureSums[index] += vals[i];
				weightedFeatureSum += val;
			}
			if(debug) {
				log.debug("Adding feature "+index+" val: "+val);
			}
		}
		Assert.a(!Double.isNaN(ret));
		return ret;
	}
}