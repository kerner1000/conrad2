package calhoun.analysis.crf.io;

import java.io.IOException;
import java.io.Serializable;
import java.util.Iterator;
import java.util.List;


/** interface to classes that handle reading and writing of input data for the CRF.  The 'location' of the input
data is passed in via a text string.  The interpretation of this text string is left up to the 
<code>InputHandler</code> implementation.  Commonly it will be a directory path, a file path, or some sort of 
configuration string.  The <code>InputHandler</code> returns {@link InputSequence}s and {@link TrainingSequence}s
to the engine.<p>
The <code>InputHandler</code> is also reponsible for writing out data.  This is necessary for subsetting and other
partitioning utilities.  As with reading, an implementation dependent location string is used to specify where the
data should be written.  When writing data, it is safe for the input handler to assume that the sequences it receives
for writing are in the same format as those it created during reading.
*/
public interface InputHandler extends Serializable {
	
	/** returns the training data read from the specified location.  Training data includes input data and 
	hidden sequences.  The result is returned as a <code>Iterator</code> so algorithms are not forced to hold
	all of the training data at once (although most will).  The interpretation of
	the location string is dependent on the particular <code>InputHandler</code> implementation used.
	@param location string location of the data.  Meaning is implementation dependent.
	@return a list of training sequences
	@throws IOException if there is a problem reading the data
	*/
	List<? extends TrainingSequence<?>> readTrainingData(String location, boolean predict) throws IOException;
	List<? extends TrainingSequence<?>> readTrainingData(String location) throws IOException;
	/** returns the input data read from the specified location.  The result is returned as an <code>Iterator</code> 
	because the inference algorithms can predict on the sequences one at a time. The interpretation of
	the location string is dependent on the particular <code>InputHandler</code> implementation used.
	@param location string location of the data.  Meaning is implementation dependent.
	@return an iterator over input sequences
	@throws IOException if there is a problem reading the data
	*/
	Iterator<? extends InputSequence<?>> readInputData(String location) throws IOException;

	/** writes training data to the specified location.  Training data includes input data and 
	hidden sequences.  The interpretation of
	the location string is dependent on the particular <code>InputHandler</code> implementation used.
	@param location string location of the data.  Meaning is implementation dependent.
	@param data a list of training sequences to write out.
	@throws IOException if there is a problem reading the data
	*/
	void writeTrainingData(String location, List<? extends TrainingSequence<?>> data) throws IOException;

	/** writes input data to the specified location.  The interpretation of
	the location string is dependent on the particular <code>InputHandler</code> implementation used.
	@param location string location of the data.  Meaning is implementation dependent.
	@param data an iterator over input sequences
	@throws IOException if there is a problem reading the data
	*/
	void writeInputData(String location, Iterator<? extends InputSequence<?>> data) throws IOException;
}
