package calhoun.analysis.crf.io;

import java.util.Collection;

import calhoun.util.Assert;

/** an input sequence where the elements of the sequence are boolean values. */
public class InputSequenceBoolean implements InputSequence<Boolean> {
	boolean[] data;
	
	/** default constructor */
	public InputSequenceBoolean() { }
	
	/** creates an input sequence using data from this boolean array 
	 * @param a a boolean array containing the input sequence */
	public InputSequenceBoolean(boolean[] a) {
		data = a;
	}
	
	/** gets the boolean array containing the data for this sequence.
	 * @return the boolean array containing the data for this sequence
	 */
	public boolean[] getData() {
		return data;
	}
	
	public Boolean getX(int ix) {
		return data[ix];
	}

	public int length() {
		return data == null ? 0 : data.length;
	}

	public InputSequence<?> getComponent(String name) {
		throw new UnsupportedOperationException();
	}

	public Collection<String> listComponents() {
		throw new UnsupportedOperationException();
	}

	public InputSequence<Boolean> subSequence(int start, int end) {
		Assert.a(start >= 1);
		Assert.a(end <= this.length());
		Assert.a(start <= end);
		boolean[] newdata = new boolean[end-start+1];
		for (int j=0; j<(end-start+1); j++) {
			newdata[j] = data[j+start-1];
		}
		InputSequenceBoolean B = new InputSequenceBoolean(newdata);
		return B;
	}
}
