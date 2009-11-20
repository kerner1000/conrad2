package calhoun.analysis.crf.io;

import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Arrays;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

public class AllIntergenicHiddenStateReader implements TrainingSequenceIO {
	private static final long serialVersionUID = 4413724139445660884L;
	private static final Log log = LogFactory.getLog(GTFInputInterval13.class);
	
	public void readTrainingSequences(Object location, List<TrainingSequence<Map<String, Object>>> seqs) throws IOException {
		for(TrainingSequence<Map<String, Object>> seq : seqs) {
			int[] states = new int[seq.length()];
			Arrays.fill(states, 0);
			seq.setY(states);
		}
		
	}
	public void writeTrainingSequences(Object location, Iterator<int[]> data) throws IOException {
		
	}
}
