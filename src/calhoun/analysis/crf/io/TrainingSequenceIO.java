package calhoun.analysis.crf.io;

import java.io.IOException;
import java.io.Serializable;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

/** Interface for reading a training sequence from the specified location. */
public interface TrainingSequenceIO extends Serializable {
	/** reads training sequences from the specified location. 
	 * @param location location of the training sequences.  Meaning is implementation dependent.
	 * @param seqs a list of the input sequences to which hidden sequences should be added.  Each hidden sequence is an array of integers containing hidden state indices.  */
	void readTrainingSequences(Object location, List<TrainingSequence<Map<String, Object>>> seqs) throws IOException;  

	/** writes training sequences to the specified location.
	 * 
	 * @param location location of the training sequences.  Meaning is implementation dependent.
	 * @param data iterator over the hidden states to write.
	 */
	void writeTrainingSequences(Object location, Iterator<int[]> data) throws IOException;  
}
