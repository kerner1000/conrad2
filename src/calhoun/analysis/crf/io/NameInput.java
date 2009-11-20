package calhoun.analysis.crf.io;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Writer;
import java.util.Map;

/** an input component that reads in a name for each sequence.  Creates an {@link NameInputSequence}.
 * Can be used as a regular input component or as part of an interleaved input.
 */
public class NameInput extends InterleavedInputComponentBase {
	private static final long serialVersionUID = -2944973705715162476L;
	
	public boolean read(BufferedReader r, Map<String, InputSequence<?>> output) throws IOException {
		String data = r.readLine();
		if(data == null)
			return false;
		output.put(name, new NameInputSequence(data));
		return true;
	}

	public void write(Writer w, Map<String, ? extends InputSequence<?>> data) throws IOException {
		w.write(((NameInputSequence) data.get(name)).getX(0));
		w.write('\n');
	}
}
