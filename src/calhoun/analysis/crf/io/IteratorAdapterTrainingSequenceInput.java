package calhoun.analysis.crf.io;

import java.util.Iterator;

/** utility class to convert an iterator over training sequences into an iterator over the underlying 
 * input sequences
 * @param <T> the type of the underlying input sequence
 */
public class IteratorAdapterTrainingSequenceInput<T> implements Iterator<InputSequence<? extends T>> {

	Iterator<? extends TrainingSequence<T>> trainingIterator;
	
	/** constructs a new iterator which will extract the input sequence from the training sequences
	 * @param trainingIterator the iterator over training sequences to extract the input sequence from
	 */
	public IteratorAdapterTrainingSequenceInput(Iterator<? extends TrainingSequence<T>> trainingIterator) {
		this.trainingIterator = trainingIterator;		
	}
	
	public boolean hasNext() {
		return trainingIterator.hasNext();
	}

	public InputSequence<? extends T> next() {
		TrainingSequence<T> trainingSeq = trainingIterator.next();
		return trainingSeq.getInputSequence();
	}

	public void remove() {
		throw new UnsupportedOperationException();
	}
}
