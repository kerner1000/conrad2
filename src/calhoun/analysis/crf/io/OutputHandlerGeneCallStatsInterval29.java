package calhoun.analysis.crf.io;

import java.io.IOException;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.util.Assert;

public class OutputHandlerGeneCallStatsInterval29 extends OutputHandlerGeneCallStats {
	private static final long serialVersionUID = 5955650475684210702L;
	private static final Log log = LogFactory.getLog(OutputHandlerGeneCallStatsInterval29.class);
	
	public void writeTestOutput(InputSequence<?> sequence, int[] truePath, int[] hiddenStates) throws IOException {
		truePath = SequenceConverter.convertSeqFromInterval29ToInterval13(truePath);
		hiddenStates = SequenceConverter.convertSeqFromInterval29ToInterval13(hiddenStates);
		//log.info("Checking hidden state validity");
		for (int i = 0; i < hiddenStates.length; i++) {
			Assert.a(hiddenStates[i] <= 12, "hiddenState[" + i + "] is " + hiddenStates[i]);
		}
		super.writeTestOutput(sequence, truePath, hiddenStates);
	}
}
