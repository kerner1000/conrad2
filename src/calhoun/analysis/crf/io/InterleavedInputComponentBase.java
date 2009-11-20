package calhoun.analysis.crf.io;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import calhoun.util.Assert;

/** base class for <code>InterleavedInputComponent</code>s.  Implements the read and write methods
 * from InputComponentIO by calling the Reader and Writer based methods of InterleavedInputComponent.
 * This allows an InterleavedInputComponent to implement just the reader and writer methods but be usable
 * in any input handler.
 * <p>
 * Also implements a name field.  The name is the name of the input sequence which will be read and written by
 * the input component.  This is the name that will identify the returned input sequence as part of a 
 * composite input.  The name defaults to "default".
 */
public abstract class InterleavedInputComponentBase implements InterleavedInputComponent {
	private static final long serialVersionUID = 2972885935425621520L;

	String name = "default";

	/** returns the name of this component
	 *  @return name of the input sequence read by this component */
	public List<String> getComponentNames() {
		return Collections.singletonList(name);
	}

	/** the name of the input sequence read and written by this component
	 * @return returns the name for the input sequence this component reads and writes
	 */
	public String getName() {
		return name;
	}

	/** sets the name of the input sequence read and written by this component
	 * @param name the name for the input sequence this component reads and writes
	 */
	public void setName(String name) {
		this.name = name;
	}

	public void readInputSequences(String location, List<Map<String, InputSequence<?>>> inputs) throws IOException {
		BufferedReader r = new BufferedReader(new FileReader(new File(location)));

		try {
			Map<String, InputSequence<?>> seq = new HashMap<String, InputSequence<?>>();
			if(inputs.size() ==  0) {
				// If the array of inputs is empty, fill it using the reader
				while(true) {
					seq = new HashMap<String, InputSequence<?>>();
					boolean success = read(r, seq);
					if(success)
						inputs.add(seq);
					else
						break;
				}
			}
			else {
				// If inputs already exist, add this input to them.  Ensure the count remains the same.
				for(Map<String, InputSequence<?>> input : inputs) {
					boolean success = read(r, input);
					Assert.a(success);
				}
				Assert.a(!read(r, seq));
			}
		}
		finally {
			r.close();
		}
	}

	public void writeInputSequences(String location, List<? extends Map<String, ? extends InputSequence<?>>> inputComponents) throws IOException {
		BufferedWriter w = new BufferedWriter(new FileWriter(new File(location)));
		try {
			for(Map<String, ? extends InputSequence<?>> seq : inputComponents) {
				write(w, seq);
			}
		}
		finally {
			w.close();
		}
	}
}
