package calhoun.analysis.crf.io;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Reader;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import calhoun.util.Assert;

/** an {@link InputHandler} for handling input files that consist of multiple different sequences
 * interleaved together in a file.  This input handler has a list of {@link InterleavedInputComponent}s.
 * When reading in a file, this InputHandler opens a reader on the file and passes the reader to each
 * {@link InterleavedInputComponent} in turn for each sequence.  Training data is assumed to be the first 
 * line of each sequence, using an {@link IntInput} to encode the hidden states.
 * <p>
 * This input handler is useful for test data that contains multiple inputs in a file along with the training data.
 * It is included to support backwards compatibility with the old input format.
 * <p>
 * This input handler can also work with "literal" input, where the location string that is passed in is not a file
 * name, but the actual input data.  This is used frequently to pass small volumes of data in unit tests.
 */
public class InputHandlerInterleaved implements InputHandler {
	private static final long serialVersionUID = -2969140424776995686L;
	
	List<InterleavedInputComponent> components;
	IntInput hiddenStateReader = new IntInput();
	boolean locationIsLiteral = false;
	boolean singleComponent = false;
	String componentName;

	/** creates a new input handler, usually to be configured from an XML file
	 */
	public InputHandlerInterleaved() {
	}
	
	/** creates a new input handler, containing a single {@link InterleavedInputComponent}
	 * @param base the single input component which is contained in the input file
	 */
	public InputHandlerInterleaved(InterleavedInputComponent base) {
		this.singleComponent = true;
		components = Collections.singletonList(base);
	}
	
	/** creates a new input handler, containing a single {@link InterleavedInputComponent}
	 * @param base the single input component which is contained in the input file
	 * @param locationIsLiteral if true then the location string passed in to the read commands
	 * is the actual input data.  Otherwise, it is the location of a file from which to read the data.
	 */
	public InputHandlerInterleaved(InterleavedInputComponent base, boolean locationIsLiteral) {
		this(base);
		this.locationIsLiteral = locationIsLiteral;
	}
	
	public Iterator<? extends InputSequence<?>> readInputData(String location) throws IOException {
		throw new UnsupportedOperationException();
	}
	
	public List<? extends TrainingSequence<?>> readTrainingData(String location) throws IOException {
		return readTrainingData(location, false);
	}
	
	public List<? extends TrainingSequence<?>> readTrainingData(String location, boolean predict) throws IOException {
		// XXX: Predict semantics are unusual here
		Reader reader = locationIsLiteral ? new StringReader(location) : new FileReader(new File(location)); 
		BufferedReader r = new BufferedReader(reader);
		List<TrainingSequence<?>> ret = new ArrayList<TrainingSequence<?>>();
		try {
			while(r.ready()) {
				int[] data = hiddenStateReader.readSequence(r);
				if(data == null) {
					break;
				}
				InputSequence<?> inputSeq = null;
				
				Map<String, InputSequence<?>> seq = new HashMap<String, InputSequence<?>>();
				for(InterleavedInputComponent comp : components) {
					boolean success = comp.read(r, seq);
					Assert.a(success == true, "Not all components of a composite input sequence were present.");
				}
				if(singleComponent) {
					Assert.a(seq.size() == 1);
					Map.Entry<String, InputSequence<?>> entry = seq.entrySet().iterator().next();
					componentName = entry.getKey();
					inputSeq = entry.getValue();
				}
				else {
					inputSeq = new InputSequenceComposite(seq);
				}
				ret.add(new TrainingSequence(inputSeq, data));
			}
		}
		finally {
			r.close();
		}
		return ret;
	}

	public void writeInputData(String location, Iterator<? extends InputSequence<?>> data) throws IOException {
		throw new UnsupportedOperationException();
	}

	public void writeTrainingData(String location, List<? extends TrainingSequence<?>> data) throws IOException {
		BufferedWriter w = new BufferedWriter(new FileWriter(new File(location)));
		try {
			for(TrainingSequence<?> seq : data) {
				hiddenStateReader.writeSequence(w, seq.getY());
				if(singleComponent) {
					Map<String, InputSequence<?>> componentSeqs = new HashMap<String, InputSequence<?>>();
					componentSeqs.put(componentName, seq.getInputSequence());
					for(InterleavedInputComponent comp : components) {
						comp.write(w, componentSeqs);
					}
				}
				else {
					InputSequenceComposite compSeq = (InputSequenceComposite) seq.getInputSequence();
					Map<String, InputSequence<?>> componentSeqs = compSeq.getComponents();
					for(InterleavedInputComponent comp : components) {
						comp.write(w, componentSeqs);
					}
				}
			}
		}
		finally {
			w.close();
		}
	}

	/** gets the current set of input components configured for this input handler.
	 * @return returns the interleaved input components that make up the file.
	 */
	public List<InterleavedInputComponent> getComponents() {
		return components;
	}

	/** sets the current set of input components configured for this input handler.
	 * @param components sets the interleaved input components that make up the file.
	 */
	public void setComponents(List<InterleavedInputComponent> components) {
		this.components = components;
	}

	/** gets the meaning of the input location string.
	 * @return true to indicate whether the input data will come in as a file
	 * or through the location string.
	 */
	public boolean isLocationIsLiteral() {
		return locationIsLiteral;
	}

	/** sets the meaning of the input location string.
	 * @param literal set locationIsLiteral to indicate whether the input data will come in as a file
	 * or through the location string.
	 */
	public void setLocationIsLiteral(boolean literal) {
		this.locationIsLiteral = literal;
	}
	
	/** gets the reader used to read in results for training data.
	 * @return the {@link TrainingSequenceIO} used to read in the hidden sequences for training
	 */
	public IntInput getHiddenStateReader() {
		return hiddenStateReader;
	}

	/** sets the reader used to get hidden sequences.  Must be set to read in training data.
	 * @param hiddenStateReader the reader that will be used to access hidden states
	 */
	public void setHiddenStateReader(IntInput hiddenStateReader) {
		this.hiddenStateReader = hiddenStateReader;
	}	
}
