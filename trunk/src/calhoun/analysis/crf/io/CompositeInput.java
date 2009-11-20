package calhoun.analysis.crf.io;

import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import calhoun.analysis.crf.Conrad;

/** a legacy input used in the original XML configuration.  No longer used. <p>
 * Represents a composite input made of several different sequences.  This is useful when putting together different
 * features that may take different inputs. 
 * If the input components are created without file entries, then the input is read from a single file and all of the components are included in the single file.
 * If the input components are created with file names, then each input component is read from its own file. 
 */
@Deprecated
public class CompositeInput implements Serializable {
	private static final long serialVersionUID = -5179361895471860850L;

	List<InputComponent> components;
	String hiddenSequenceFile = "hidden.dat";
	TrainingSequenceIO hiddenSequenceReader = new IntInput();
	
	public static class InputComponent implements Serializable {
		private static final long serialVersionUID = -652343588509377034L;
		String name;
		InterleavedInputComponentBase inputSequence;
		String filename;
		
		public String getFilename() {
			return filename;
		}
		public void setFilename(String filename) {
			this.filename = filename;
		}
		public InterleavedInputComponentBase getInputSequence() {
			return inputSequence;
		}
		public void setInputSequence(InterleavedInputComponentBase inputSequence) {
			this.inputSequence = inputSequence;
		}
		public String getName() {
			return name;
		}
		public void setName(String name) {
			this.name = name;
		}
	}
	
	public void setComponents(List<InputComponent> components) {
		this.components = components;
	}

	public String getHiddenSequenceFile() {
		return hiddenSequenceFile;
	}

	public void setHiddenSequenceFile(String hiddenSequenceFile) {
		this.hiddenSequenceFile = hiddenSequenceFile;
	}

	public TrainingSequenceIO getHiddenSequenceReader() {
		return hiddenSequenceReader;
	}
	
	public void setHiddenSequenceReader(TrainingSequenceIO hiddenSequenceReader) {
		this.hiddenSequenceReader = hiddenSequenceReader;
	}
	
	/** legacy {@link InputHandler} that mirrors the original XML config file behavior.  
	 * There is a small amount of support in {@link Conrad} to set this up to be truly backwards compatible.
	*/
	public static class LegacyInputHandler implements InputHandler {
		private static final long serialVersionUID = -5031580949745955877L;

		Object inputComponent;
		InputHandler handler;
		
		/** special constructor used by {@link Conrad} to pass the old style input config to the legacy input handler.
		 * @param inputComponent the originally configured inputFormat bean.
		 */
		public LegacyInputHandler(Object inputComponent) {
			this.inputComponent = inputComponent;
		}
		
		public Iterator<? extends InputSequence<?>> readInputData(String inputLocation) throws IOException {
			initInputHandler(inputLocation);
			return handler.readInputData(inputLocation);
		}
		
		public List<? extends TrainingSequence<?>> readTrainingData(String location) throws IOException {
			return readTrainingData(location, false);
		}
		
		public List<? extends TrainingSequence<?>> readTrainingData(String inputLocation, boolean predict) throws IOException {
			initInputHandler(inputLocation);
			return handler.readTrainingData(inputLocation, predict);
		}
		
		public void writeTrainingData(String location, List data) throws IOException {
			initInputHandler(location);
			handler.writeTrainingData(location, data);
		}

		public void writeInputData(String location, Iterator data) throws IOException {
			initInputHandler(location);
			handler.writeInputData(location, data);
		}

		/** default behavor.  If a single input component is listed, always assume it is a file.  If a composite is specified,
		 * choose interleaved or directory based on the value of the inputLocation parameter.
		 */
		void initInputHandler(final String inputLocation) {
			if(handler != null)
				return;
			if(!(inputComponent instanceof CompositeInput)) {
				handler = new InputHandlerInterleaved((InterleavedInputComponent) inputComponent);
				return;
			}
			final CompositeInput compositeInput = (CompositeInput) inputComponent;
			File file = new File(inputLocation);
			if(file.isDirectory()) {
				InputHandlerDirectory d = new InputHandlerDirectory();
				d.setHiddenSequenceFile(compositeInput.getHiddenSequenceFile());
				d.setHiddenStateReader(compositeInput.getHiddenSequenceReader());

				Map<String, InputComponentIO> components = new HashMap();
				for(InputComponent comp : compositeInput.components) {
					comp.getInputSequence().setName(comp.getName());
					components.put(comp.getFilename(), comp.getInputSequence());
				}
				d.setInputReaders(components);
				handler = d;
			}
			else {
				InputHandlerInterleaved interleaved = new InputHandlerInterleaved();
				interleaved.setHiddenStateReader((IntInput)compositeInput.getHiddenSequenceReader());
				List<InterleavedInputComponent> components = new ArrayList();
				for(InputComponent comp : compositeInput.components) {
					comp.getInputSequence().setName(comp.getName());
					components.add(comp.getInputSequence());
				}
				interleaved.setComponents(components);
				handler = interleaved;
			}
		}
	}
}
