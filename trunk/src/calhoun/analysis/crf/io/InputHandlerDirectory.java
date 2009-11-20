package calhoun.analysis.crf.io;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import calhoun.util.Util;

/** an {@link InputHandler} used when the input is in several files within a single directory.  A single {@link InputComponentIO} is used for each
 * file.  A map associates each file name with its {@link InputComponentIO}.  For training, hidden sequences are stored in a separate file in the directory whose name is set with the hiddenSequenceFile property.
 * For this {@link InputHandler}, the location passed is the path to the directory containing the input data.
  */
public class InputHandlerDirectory extends InputHandlerBase {
	private static final long serialVersionUID = -2969140424776995686L;
	
	Map<String, InputComponentIO> inputReaders;
	TrainingSequenceIO hiddenStateReader;
	String hiddenSequenceFile = "hidden.dat";

	public Iterator<? extends InputSequence<?>> readInputData(String location) throws IOException {
		List<Map<String, InputSequence<?>>> inputs = new ArrayList();

		// Read in all of the inputs
		for(Map.Entry<String, InputComponentIO> entry : inputReaders.entrySet()) {
			entry.getValue().readInputSequences(new File(location, entry.getKey()).getPath(), inputs);
		}

		return createCompositeInput(inputs);
	}

	public List<? extends TrainingSequence<?>> readTrainingData(String location) throws IOException {
		return readTrainingData(location, false);
	}
	
	public List<? extends TrainingSequence<?>> readTrainingData(String location, boolean predict) throws IOException {
		String trainingLocation = new File(location, hiddenSequenceFile).getPath();

		return readTrainingData(location, trainingLocation, hiddenStateReader, predict);
	}

	public void writeInputData(String location, Iterator<? extends InputSequence<?>> data) throws IOException {
		// Collect all the values from the iterator into a list
		// Then for each composite, separate it into a map of its component pieces for handing to the IO class
		List<Map<String, InputSequence<?>>> compList = new ArrayList<Map<String, InputSequence<?>>>();
		Util.addAll(compList, new IteratorAdapterInputComponent(data));

		for(Map.Entry<String, InputComponentIO> entry : inputReaders.entrySet()) {
			entry.getValue().writeInputSequences(new File(location, entry.getKey()).getPath(), compList);
		}
	}

	public void writeTrainingData(String location, List<? extends TrainingSequence<?>> data) throws IOException {
		writeInputData(location, new IteratorAdapterTrainingSequenceInput(data.iterator()));

		List<int[]> trainingSeqs = new ArrayList<int[]>();
		for(TrainingSequence<?> t : data) {
			trainingSeqs.add(t.getY());
		}
		
		hiddenStateReader.writeTrainingSequences(new File(location, hiddenSequenceFile).getPath(), trainingSeqs.iterator());
	}

	/** gets the reader used to read in results for training data.
	 * @return the {@link TrainingSequenceIO} used to read in the hidden sequences for training
	 */
	public TrainingSequenceIO getHiddenStateReader() {
		return hiddenStateReader;
	}

	/** sets the reader used to get hidden sequences.  Must be set to read in training data.
	 * @param hiddenStateReader the reader that will be used to access hidden states
	 */
	public void setHiddenStateReader(TrainingSequenceIO hiddenStateReader) {
		this.hiddenStateReader = hiddenStateReader;
	}
	
	/** gets the readers used to read in input sequences.  Must be set before any of the <code>read</code> methods are called.
	 * @return the reader used to read in input sequences.
	 */
	public Map<String, InputComponentIO> getInputReaders() {
		return inputReaders;
	}
	
	/** sets the readers used to read in input sequences.  Must be set before any of the <code>read</code> methods are called.
	 * the value is a map that associates filenames within the directory to input components.
	 * @param inputReader the reader used to read in input sequences.
	 */
	public void setInputReaders(Map<String, InputComponentIO> inputReader) {
		this.inputReaders = inputReader;
	}
	
	/** gets the name of the hidden sequence file.  This is the name of the file within the directory where training data will be located.
	 * @return the name of the hidden sequence file.
	 */
	public String getHiddenSequenceFile() {
		return hiddenSequenceFile;
	}

	/** sets the name of the hidden sequence file.  This is the name of the file within the directory where training data will be located.
	 * @param hiddenSequenceFile the name of the hidden sequence file within the input directory.
	 */
	public void setHiddenSequenceFile(String hiddenSequenceFile) {
		this.hiddenSequenceFile = hiddenSequenceFile;
	}
}
