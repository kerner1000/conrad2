/**
 * 
 */
package calhoun.analysis.crf.solver;

import java.util.Iterator;

/** Utility class that implements a FIFO buffer with constant time adds, removes,
 * and accesses to any position using a circular buffer in an array.
 */
public final class RecyclingBuffer<T> {
	public T[] array;
	public int length;
	public int currentStart = length;

	public RecyclingBuffer(final T[] arg) {
		this.array = arg;
		length = array.length;
	}
	
	public final T get(final int i) {
		if(i >= length)
			throw new ArrayIndexOutOfBoundsException(i);
		return array[(currentStart+i)%length];
	}

	public final void set(final int i, final T val) {
		if(i >= length)
			throw new ArrayIndexOutOfBoundsException(i);
		array[(currentStart+i)%length] = val;
	}

	public final T addFirst(final T t) {
		currentStart = (currentStart+length-1)%length; 
		T ret = array[currentStart];
		array[currentStart] = t;
		return ret;
	}
	
	public final int getTotalSize() {
		return length;
	}
	
	final class RecyclingBufferIterator implements Iterator<T> {
		int current = 0;
		public final boolean hasNext() {
			return current < length;
		}
		public final T next() {
			T ret = get(current);
			++current;
			return ret;
		}

		public final void remove() {
			throw new UnsupportedOperationException();
		}
	}
}