package calhoun.analysis.crf.io;

import java.io.BufferedReader;
import java.io.IOException;

public class IntInputInterval29 extends IntInput {
	private static final long serialVersionUID = -2301269002866688631L;

	public int[] readSequence(BufferedReader r) throws IOException {
		int[] ret;
		ret = super.readSequence(r);
		ret = SequenceConverter.convertSeqFromInterval13ToInterval29(ret);
		return ret;
	}
}
