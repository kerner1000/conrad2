package calhoun.analysis.crf.solver.check;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.FeatureList;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.util.Assert;
import calhoun.util.ErrorException;
import cern.colt.matrix.DoubleMatrix1D;

/** Implementation of the FeatureList interface.  An interface is used as the public interface so that
 * features cannot get access to internal values.  This class should be used by solvers to access the
 * feature values that are returned. */
public class ArrayFeatureList implements FeatureList {
	private static final Log log = LogFactory.getLog(ArrayFeatureList.class);
	private static final boolean debug= log.isDebugEnabled();

	// How much to grow the arrays if necessary.  Set low because it shouldn't usually happen.
	static final float GROW_RATE = 1.5f;
	
	private ModelManager manager;
	public boolean valid;
	public int[] indices;
	public double[] values;
	public int currentSize;
	
	/** Creates arrays of the given size to hold the features.  Arrays will grow automatically if required. */
	public ArrayFeatureList(ModelManager manager) {
		this.manager = manager;
		int length = manager.getNumFeatures()+1;
		if(debug) {
			log.debug("Creating feature list of size "+length);
		}
		indices = new int[length];
		values = new double[length];
		clear();
	}
	
	/** Adds the feature onto the array at the next postion.  Expands the arrays if necessary
	 */
	public void addFeature(int index, double value) {
		if(Double.isInfinite(value) || Double.isNaN(value)) {
			Assert.a(false, "Feature #"+index+" returned "+value);
		}
		//if(debug) {
		//	log.debug(String.format("Add feat %d = %.2f", index, value));
		//}
		/*if(value == 0.0) {
			return;
		}*/
		if(currentSize == values.length) {
			int oldLength = values.length;
			int newLength = (int) (GROW_RATE * oldLength + 1);
			log.warn("Expanding feature list from "+currentSize+" to "+newLength);
			double[] newValues = new double[newLength];
			int[] newIndices = new int[newLength];
			System.arraycopy(values, 0, newValues, 0, oldLength);
			System.arraycopy(indices, 0, newIndices, 0, oldLength);
			values = newValues;
			indices = newIndices;
		}
		indices[currentSize] = index;
		values[currentSize] = value;
		currentSize++;
	}
	
	/** Returns the invalid flag. */
	public boolean isValid() {
		return valid;
	}

	/** Invalidates results. */
	public void invalidate() {
		valid = false;
	}

	public void evaluateNode(InputSequence seq, int pos, int state) {
		clear();
		manager.evaluateNode(seq, pos, state, this);
	}

	public void evaluateEdge(InputSequence seq, int pos, int prevState, int state) {
		clear();
		manager.evaluateEdge(seq, pos, prevState, state, this);
	}

	public void evaluateNodeLength(InputSequence seq, int pos, int len, int state) {
		clear();
		manager.evaluateNodeLength(seq, pos, len, state, this);
	}

	public void evaluateEdgeLength(InputSequence seq, int pos, int len, int prevState, int state) {
		clear();
		manager.evaluateEdgeLength(seq, pos, len, prevState, state, this);
	}

	/** Resets the structure */
	public void clear() {
		valid = true;
		currentSize = 0;
	}

	public int getIndex(int index) {
		Assert.a(index < currentSize);
		return indices[index];
	}
	
	public double getValue(int index) {
		Assert.a(index < currentSize);
		return values[index];
	}
	
	/** Returns the internal array of feature ids that were added.  Must use the size() function to find out how many array eleents are valid. */
	public int[] getIndices() {
		return indices;
	}
	
	/** Returns the internal array of feature values that were added.  Must use the size() function to find out how many array eleents are valid. */
	public double[] getValues() {
		return values;
	}

	/** Returns the number of features added since the last clear */
	public int size() {
		return currentSize;
	}

	public void updateExpectations(DoubleMatrix1D expects, double multiplier) {
		int size = size();
		for(int feat = 0; feat < size; ++feat) {
			int index = indices[feat];
			double val = values[feat] * multiplier;
			//log.debug("Setting expect: "+val);
			expects.setQuick(index, val + expects.getQuick(index));
			if(Double.isNaN(val)) {
				throw new ErrorException("Overflow");
			}
		}
	}
}
