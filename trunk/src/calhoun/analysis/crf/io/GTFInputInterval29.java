package calhoun.analysis.crf.io;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.io.GTFInputInterval13.CDS;
import calhoun.util.Assert;
import calhoun.util.FileUtil;

public class GTFInputInterval29 extends GTFInputInterval13 {
	private static final long serialVersionUID = 4413724139445660884L;
	private static final Log log = LogFactory.getLog(GTFInputInterval29.class);
	
	public void readTrainingSequences(Object location, List<TrainingSequence<Map<String, Object>>> seqs) throws IOException {
		super.readTrainingSequences(location, seqs);
		
		for (Iterator iter = seqs.iterator(); iter.hasNext();) {
			TrainingSequence<Map<String, Object>> seq = (TrainingSequence<Map<String, Object>>)iter.next();
			int[] states = seq.getY();
			states = SequenceConverter.convertSeqFromInterval13ToInterval29(states);
			seq.setY(states);
		}
	}
}
