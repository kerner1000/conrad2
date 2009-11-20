package calhoun.analysis.crf.io;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Writer;
import java.util.List;
import java.util.Map;

import calhoun.util.Assert;

/** reads in an input consisting of a string.  Can be used as a standalone
 * input component or part of an interleaved input.
 */
public class StringInput extends InterleavedInputComponentBase  {
	private static final long serialVersionUID = -2944973705715162476L;

	public boolean read(BufferedReader r, Map<String, InputSequence<?>> output) throws IOException {
		String data = r.readLine();
		if(data == null)
			return false;
		output.put(name, new InputSequenceCharacter(data));
		return true;
	}

	public void write(Writer w, Map<String, ? extends InputSequence<?>> data) throws IOException {
		InputSequenceCharacter comp = (InputSequenceCharacter) data.get(name);
		Assert.a(comp != null, "No component found in input: "+name+".  Entry is "+data.keySet().iterator().next());
		w.write(comp.getString());
		w.write('\n');
	}

	/** Convenience function for creating training sequences in test data. */
	public static List<? extends TrainingSequence<Character>> prepareData(String str) {
		try {
			return (List<? extends TrainingSequence<Character>>) new InputHandlerInterleaved(new StringInput(), true).readTrainingData(str) ;
		}
		catch(Exception ex) {
			throw new RuntimeException(ex);
		}
	}
}
