package calhoun.analysis.crf.io;

import java.util.Collection;

import calhoun.util.Assert;

/** an input sequence where the elements of the sequence are objects. */
public class InputSequenceObject<T> implements InputSequence<T> {
	T[] t;
	
	/** default constructor */
	public InputSequenceObject() {
	}
	
	/** creates an input sequence using data from this Object array 
	 * @param t a an object array containing the input sequence */
	public InputSequenceObject(T[] t) {
		this.t = t;
	}
	
	public T getX(int ix) {
		return t[ix];
	}

	public int length() {
		return t.length;
	}

	public InputSequence<?> getComponent(String name) {
		throw new UnsupportedOperationException();
	}

	public Collection<String> listComponents() {
		throw new UnsupportedOperationException();
	}

	public InputSequence<T> subSequence(int start, int end) {
		Assert.a(start >= 1);
		Assert.a(end <= this.length());
		Assert.a(start <= end);

		int length = end - start + 1;
		T[] ret = (T[]) new Object[length];
		for(int i= 0; i<length; ++i) {
			ret[i] = t[i+start-1];
		}
		return new InputSequenceObject<T>(ret);
	}
}
