package calhoun.analysis.crf.io;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import calhoun.util.Assert;
import calhoun.util.ErrorException;

/** reads in an input consisting of a list of ints that correspond to binary values.  Can be used as a standalone
 * input component or part of an interleaved input.
 */
public class IntInput extends InterleavedInputComponentBase implements TrainingSequenceIO {
	private static final long serialVersionUID = 4413724139445660883L;
	
	public boolean read(BufferedReader r, Map<String, InputSequence<?>> output) throws IOException {
		int[] data = readSequence(r);
		if(data == null) {
			return false;
		}
		output.put(name, new InputSequenceInt(data));
		return true;
	}
	
	public int[] readSequence(BufferedReader r) throws IOException {
		String str = r.readLine();
		if(str == null) {
			return null;
		}
		int[] data = new int[str.length()];
		try {
			for (int i = 0; i < str.length(); ++i) {
				int temp = str.charAt(i) - '0';
				if ( (temp<0) || (temp>9)) {
					temp = str.charAt(i) - 'A' + 10;
					
					if ( (temp<10) || (temp>35)) {
						temp = str.charAt(i) - 'a' + 36;
						Assert.a( (temp>=36) && (temp<62), "Offending character was '" + str.charAt(i));
					}
				}
				data[i] = temp;
			}
		} catch (NumberFormatException ex) {
			throw new ErrorException(ex);
		}
		return data;
	}
	
	public void write(Writer w, Map<String, ? extends InputSequence<?>> data) throws IOException {
		writeSequence(w, ((InputSequenceInt) data.get(name)).getData());
	}

	public void writeSequence(Writer w, int[] data) throws IOException {
		for(int i: data) {
			if (i<10) {
				w.write('0'+i);
			} else if (i<36) {
				w.write('A'+(i-10));
			} else if (i<62) {
				w.write('a'+(i-36));
			} else { throw new IOException(); }
		}
		w.write('\n');
	}

	public void readTrainingSequences(Object location, List<TrainingSequence<Map<String, Object>>> seqs) throws IOException {
		BufferedReader r = new BufferedReader(new FileReader(new File((String) location)));
		Iterator<TrainingSequence<Map<String, Object>>> it = seqs.iterator();
		while(r.ready()) {
			int[] data = readSequence(r);
			if(data == null)
				break;
			TrainingSequence<Map<String, Object>> seq = it.next();
			seq.setY(data);
		}
	}

	public void writeTrainingSequences(Object location, Iterator<int[]> data) throws IOException {
		BufferedWriter w = new BufferedWriter(new FileWriter(new File((String) location)));
		while(data.hasNext()) {
			int[] seq = data.next();
			writeSequence(w, seq);
		}
		w.close();
	}

	/** Convenience function for creating training sequences in test data. */
	public static List<? extends TrainingSequence<?>> prepareData(String str) throws Exception {
		return new InputHandlerInterleaved(new IntInput(), true).readTrainingData(str) ;
	}
}
