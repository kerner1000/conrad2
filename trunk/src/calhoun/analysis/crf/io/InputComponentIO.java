package calhoun.analysis.crf.io;

import java.io.IOException;
import java.io.Serializable;
import java.util.List;
import java.util.Map;

/** interface to read or write {@link InputSequence} components to or from a file.  A single <code>InputComponentIO</code> is responsible for 
 * reading and writing data at a single location, but may work with multiple components of the input sequence.
 */
public interface InputComponentIO extends Serializable {
	
	/** A list of names of the components of the InputSequence created by this reader. 
	 * @return a list of input sequence component names. */
	List<String> getComponentNames();
	
	/** reads input sequences from this location.  In most cases the location will be a file and this object will read in one or more components of the
	 * overall input from that file.  Alist of inputs is returned.  Each input consists of a set of key-value pairs, where the key is the name of the 
	 * input component and the value is an (@link InputSequence) object.  If this is not the first component to be loaded, this list may already contain
	 * entries for the input.  
	 * @param location the location of the data to read.  Meaning of the location is implementation dependent, but will usually be a file name.
	 * @param inputs a list of input sequences.  For each input sequence, a map is returned that maps component names to their
	 * associated {@link InputSequence} objects.
	 */
	void readInputSequences(String location, List<Map<String, InputSequence<?>>> inputs) throws IOException;  

	/** writes input sequences to this location.  <code>InputComponentIO</code> optionally can implement this function to provide the ability to write input
	 * data as well as read it.  This is used by many of the data manipulation tools, such as those for subsetting and creating cross-validation sets.
	 * @param location the location of the data to write.  Meaning of the location is implementation dependent, but will usually be a file name.
	 * @param inputComponents an iterator over input sequences.  For each input sequence, a map is returned that maps component names to their
	 * associated {@link InputSequence} objects.  All components of the input are passed, and this object is reponsible for knowing which components
	 * it should be writing.
	 */
	void writeInputSequences(String location, List<? extends Map<String, ? extends InputSequence<?>>> inputComponents)  throws IOException;  
}
