package calhoun.analysis.crf.io;

import java.util.Iterator;
import java.util.Map;

/** utility adapter class for creating an iterator of input components from an {@link InputSequenceComposite} iterator.
 *  */
public class IteratorAdapterInputComponent implements Iterator<Map<String, InputSequence<?>>> {
	
	Iterator<? extends InputSequence<?>> parent;

	/** constructs an iterator that separates out the input components of the parent iterator
	 * @param parent the iterator over composite input sequences
	 */
	public IteratorAdapterInputComponent(Iterator<? extends InputSequence<?>> parent) {
		this.parent = parent;
	}
	
	public boolean hasNext() {
		return parent.hasNext();
	}

	public Map<String, InputSequence<?>> next() {
		InputSequenceComposite seq = (InputSequenceComposite) parent.next();
		return seq.getComponents();
	}

	public void remove() {
		throw new UnsupportedOperationException();
	}
}
