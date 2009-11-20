package calhoun.analysis.crf.io;

import java.io.IOException;
import java.io.Serializable;

/** handles outputting of results from predictions. */
public interface OutputHandler extends Serializable {

	/** sets the output location to write to.  This is called before any of the write functions are called.
	 * @param location the location to write the data to.  The meaning is implementation dependent.
	 */
	void setOutputLocation(String location);
	
	/** Writes out a set of predicted hidden states.  The location will have been specified previously with {@link #setOutputLocation}
	 * @param sequence the input sequence on which predictions were made
	 * @param hiddenStates an array of state indices for the predicted hidden states.
	 */
	void writeOutput(InputSequence<?> sequence, int[] hiddenStates) throws IOException;

	/** Writes out a set of hidden states compared with a known true path.  The location will have been specified previously with {@link #setOutputLocation}
	 * @param sequence the input sequence on which predictions were made
	 * @param truePath an array of hidden states representing the true path.
	 * @param hiddenStates an array of state indices for the predicted hidden states.
	 */
	void writeTestOutput(InputSequence<?> sequence, int[] truePath, int[] hiddenStates) throws IOException;

	/** indicates the writing of output is complete.  The handler can do any final processing and release resources.
	 */
	void outputComplete() throws IOException;
}
