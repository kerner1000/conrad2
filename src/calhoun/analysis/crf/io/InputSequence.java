package calhoun.analysis.crf.io;

import java.util.Collection;

/** represents the observed inputs.  Contains a sequence of inputs, and allows them to be retrieved by position.<p>
 * An input sequence may be made up of different components.  If it contains components the {@link #getComponent} and 
 * {@link #listComponents()} functions must be implemented.  Otherwise they can just trhow {@link UnsupportedOperationException}.
 * <p>
 * In some cases, the input may not be a sequence.  It may be a single input value that is returned for all positions.  In
 * this case, the length value can return -1 and the {@link #getX} function can return the same object at each position.
 * */
public interface InputSequence<A> {
	
	/** retrieves the input value at a position in the input sequence.
	 * @param x the index position at which to get the input.  This is a zero-based index.
	 * @return the object at this position in the input.
	 */
	A getX(int x);
	
	/** Returns the length of this sequence 
	 * @return length */
	int length() ;
	
	/**tTakes a subinterval of the input sequence with given start-end coordinates
	 * which are relative coordinates, 1-based, and inclusive.  Thus 1-10 will mean
	 * returning a new InputSequence which is the first 10 positions of the current one.<p>  
	 * An implementation that does not support subsetting should throw an {@link UnsupportedOperationException}
	 * @param start the 1-based index of the first position of the input to retrieve.
	 * @param end the 1-based index of the last position  of the input to retrieve.
	 * @return an input sequence which is a sbusequence of the original sequence.
	 *  */
	InputSequence<A> subSequence(int start, int end) ;
	
	/** For input sequences that are a composite of several different input objects, returns 
	 * a particular component of the input. For simple input sequences that just return an object, 
	 * this method should throw {@link UnsupportedOperationException}. 
	 * @param name the name of the input component to return 
	 * @return the input sequence representing this component */
	InputSequence<?> getComponent(String name);

	/** For input sequences that are a composite of several different input objects, returns a list of the names of the components in this input sequence. 
	 * For simple input sequences that just return an object, this method should throw UnsupportedOperationException. 
	 * @return a collection containing the component names */
	Collection<String> listComponents();
}
