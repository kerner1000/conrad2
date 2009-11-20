package calhoun.analysis.crf.io;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Writer;
import java.util.Map;

/** an input component that can be used with the {@link InputHandlerInterleaved}.  This is a regular
 * input component that also has input and output methods that take a BufferedReader to use when reading data
 * from the input file. */
public interface InterleavedInputComponent extends InputComponentIO {

	/** Read the contents of the sequence in from a reader.  Return false if the end of file was reached. 
	 * @param r the buffered reader from which the next sequence should be read.
	 * @param output a map to which new components for this input sequence should be added. 
	 * @return true if a sequence was read in and the map was populated.  false if the end of the sequence was reached. */
	boolean read(BufferedReader r, Map<String, InputSequence<?>> output) throws IOException;

	/** Output this sequence to the given writer.
	 * @param w a writer to which this sequence should be written
	 * @param data a map containing the data for the sequence to be written */
	void write(Writer w, Map<String, ? extends InputSequence<?>> data) throws IOException;
}
