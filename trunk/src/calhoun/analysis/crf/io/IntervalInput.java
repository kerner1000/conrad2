package calhoun.analysis.crf.io;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Writer;
import java.util.Map;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.io.IntervalInputSequence.IntervalRangeMapValue;
import calhoun.util.Assert;
import calhoun.util.RangeMap;

/** an input component used to read in values that are similar for long intervals.
 * Each input sequence is a line which looks like "<code>(164716,215910,-,1.0)(194048,218199,-,1.0)...</code>".
 * Takes this input and creates an {@link IntervalInputSequence}.
 * <p>
 * This input component is can be used in a regular input handler or as part of an {@link InputHandlerInterleaved}.
 */
public class IntervalInput extends InterleavedInputComponentBase {
	private static final long serialVersionUID = 8110670530135286052L;
	private static final Log log = LogFactory.getLog(IntervalInput.class);

	public boolean read(BufferedReader r, Map<String, InputSequence<?>> output) throws IOException {
		String temp = r.readLine();
		String inputName = temp;
		temp = r.readLine();
		int inputLength;
		try {
			inputLength = Integer.parseInt(temp);
		} catch (Exception e){
			throw new RuntimeException("Offending line was : " + temp);
		}
		log.debug("Length of " + inputName + "   is  " + inputLength);	
		
		RangeMap rmplus = new RangeMap();
		RangeMap rmminus = new RangeMap();
		
		String bigLine = r.readLine();
		// Parse lines like this: (164716,215910,-,1.0)(194048,218199,-,1.0)	
		
		String[] fields = bigLine.split("[(]");
		for (int k=1; k<fields.length; k++) {
			String[] rings = fields[k].split("[,)]");

			int start = Integer.parseInt(rings[0]);
			int stop = Integer.parseInt(rings[1]);
			char strand = rings[2].charAt(0);
			double val = Double.parseDouble(rings[3]);

			IntervalRangeMapValue irmv = new IntervalRangeMapValue(start,stop,val);
			if (strand=='+') {
				irmv.insertIntoRangeMap(rmplus);
			} else if (strand=='-') {
				irmv.insertIntoRangeMap(rmminus);			
			} else { Assert.a(false); }
		}
		InputSequence<?> inputSeq = new IntervalInputSequence(rmplus, rmminus, inputName, inputLength);
		output.put(name, inputSeq);
		return true;
	}

	public void write(Writer w, Map<String, ? extends InputSequence<?>> data) throws IOException {
		IntervalInputSequence seq = (IntervalInputSequence) data;
		w.write(seq.inputName + "\n");
		w.write("" + seq.inputLength + "\n");
		Object[] rmplusList = seq.rmplus.getObjectList().values().toArray();
		for (int i=0; i<rmplusList.length; i++) {
			w.write(((IntervalRangeMapValue) rmplusList[i]).toStringStrand("+"));
		}
		
		Object[] rmminusList = seq.rmminus.getObjectList().values().toArray();
		for (int i=0; i<rmminusList.length; i++) {
			w.write(((IntervalRangeMapValue) rmminusList[i]).toStringStrand("+"));
		}		
		w.write("\n");
	}
}
