package calhoun.analysis.crf.io;

import java.util.Collection;

import calhoun.util.Assert;

/** an input sequence where the elements of the sequence are integer values. */
public class InputSequenceInt implements InputSequence<Integer> {
	int[] data;
	
	/** default constructor */
	public InputSequenceInt() { }
	
	/** creates an input sequence using data from this int array 
	 * @param a an int array containing the input sequence */
	public InputSequenceInt(int[] a) {
		data = a;
	}
	
	/** gets the int array containing the data for this sequence.
	 * @return the int array containing the data for this sequence
	 */
	public int[] getData() {
		return data;
	}
	
	public Integer getX(int ix) {
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

	public InputSequence<Integer> subSequence(int start, int end) {
		Assert.a(start >= 1);
		Assert.a(end <= this.length());
		Assert.a(start <= end);
		int[] newdata = new int[end-start+1];
		for (int j=0; j<(end-start+1); j++) {
			newdata[j] = data[j+start-1];
		}
		InputSequenceInt B = new InputSequenceInt(newdata);
		return B;
	}
}
