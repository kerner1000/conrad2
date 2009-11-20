package calhoun.analysis.crf.io;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import calhoun.util.Util;

/** an {@link InputHandler} used when all of the input is in a single file.  A single {@link InputComponentIO} is used to read the
 * file.  For training, hidden sequences are stored in a separate file whose name is related to the input file name using a
 * {@link FilenameMapper}.  The same filename mapping is used to determine the training set file name when writing out data
 * as when reading it in.  The training file is read using a {@link TrainingSequenceIO}.<p>
 * For this {@link InputHandler}, the location passed is the path to the file containing the input data.
  */
public class InputHandlerFile extends InputHandlerBase {
	private static final long serialVersionUID = -2969140424776995686L;
	
	InputComponentIO inputReader;
	FilenameMapper mapper;
	TrainingSequenceIO hiddenStateReader;
	
	public Iterator<? extends InputSequence<?>> readInputData(String location) throws IOException {
		List<Map<String, InputSequence<?>>> inputs = new ArrayList();
		inputReader.readInputSequences(location, inputs);

		return createCompositeInput(inputs);
	}

	public List<? extends TrainingSequence<?>> readTrainingData(String location) throws IOException {
		return readTrainingData(location, false);
	}
	
	public List<? extends TrainingSequence<?>> readTrainingData(String location, boolean predict) throws IOException {
		String trainingLocation = mapper.mapFilename(new File(location)).getPath();
		return readTrainingData(location, trainingLocation, hiddenStateReader, predict);
	}

	public void writeInputData(String location, Iterator<? extends InputSequence<?>> data) throws IOException {
		// Collect all the values from the iterator into a list
		// Then for each composite, separate it into a map of its component pieces for handing to the IO class
		List<Map<String, InputSequence<?>>> compList = new ArrayList<Map<String, InputSequence<?>>>();
		Util.addAll(compList, new IteratorAdapterInputComponent(data));

		inputReader.writeInputSequences(location, compList);
	}

	public void writeTrainingData(String location, List<? extends TrainingSequence<?>> data) throws IOException {
		writeInputData(location, data.iterator());

		String trainingLocation = mapper.mapFilename(new File(location)).getPath();
		
		List<int[]> trainingSeqs = new ArrayList<int[]>();
		for(TrainingSequence<?> t : data) {
			trainingSeqs.add(t.getY());
		}
		
		hiddenStateReader.writeTrainingSequences(trainingLocation, trainingSeqs.iterator());
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
	
	/** gets the reader used to read in input sequences.  Must be set before any of the <code>read</code> methods are called.
	 * @return the reader used to read in input sequences.
	 */
	public InputComponentIO getInputReader() {
		return inputReader;
	}
	
	/** gets the reader used to read in input sequences.  Must be set before any of the <code>read</code> methods are called.
	 * @param inputReader the reader used to read in input sequences.
	 */
	public void setInputReader(InputComponentIO inputReader) {
		this.inputReader = inputReader;
	}
	
	/** the mapper used to generate the name of the hidden sequence file from the input sequence file.  
	 * Must be set to read in training data.
	 * @return the mapper used to generate the hidden sequence file name.
	 */
	public FilenameMapper getMapper() {
		return mapper;
	}

	/** the mapper used to generate the name of the hidden sequence file from the input sequence file.
	 * @param mapper the mapper used to generate the hidden sequence file name.
	 */
	public void setMapper(FilenameMapper mapper) {
		this.mapper = mapper;
	}
}
