package calhoun.util;

import java.util.Iterator;

public abstract class FilterIterator<A> implements Iterator<A> {
	Iterator<A> baseIterator;
	A current = null;
	
	
	public FilterIterator(Iterator<A> base) {
		baseIterator = base;
		advance();
	}
	
	void advance() {
		boolean valid = false;
		while(baseIterator.hasNext()) {
			current = baseIterator.next();
			if(test(current)) {
				valid = true;
				break;
			}
		}
		if(!valid)
			current = null;
	}
	
	public boolean hasNext() {
		return current != null;
	}
	
	public A next() {
		A ret = current;
		advance();
		return ret;
	}

	public void remove() {
		throw new UnsupportedOperationException();
	}
	
	public abstract boolean test(A val);
}
