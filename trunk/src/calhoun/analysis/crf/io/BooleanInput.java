package calhoun.analysis.crf.io;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Writer;
import java.util.Map;

import calhoun.util.Assert;

/** reads in an input consisting of a string of 1's and 0's that correspond to binary values.  Can be used as a standalone
 * input component or part of an interleaved input.
 */
public class BooleanInput extends InterleavedInputComponentBase {
	private static final long serialVersionUID = 4922000330136279956L;
	
	public boolean read(BufferedReader r, Map<String, InputSequence<?>> output) throws IOException {
		String str = r.readLine();
		if(str == null) {
			return false;
		}
		boolean[] data = new boolean[str.length()];
		for (int i = 0; i < str.length(); ++i) {
			char c = str.charAt(i);
			if(c == '1')
				data[i] = true;
			else {
				Assert.a(c=='0');
			}
		}
		output.put(name, new InputSequenceBoolean(data));
		return true;
	}
	
	public void write(Writer w, Map<String, ? extends InputSequence<?>> data) throws IOException {
		boolean[] boolData = ((InputSequenceBoolean) data.get(name)).getData();
		for(boolean i: boolData) {
			w.write(i ? '1':'0');
		}
		w.write('\n');
	}
}
