package calhoun.analysis.crf.test;

import calhoun.analysis.crf.Conrad;
import calhoun.util.AbstractTestCase;
import calhoun.util.FileUtil;

public class MultipleSequenceTest extends AbstractTestCase {

	public void testMultipleSequences() throws Exception {
		// Verify that we get equivalent results when evaluating a gradient for multiple sequences.
		// Train features against a set of sequences.
		Conrad crf1 = new Conrad("test/input/multipleSequence/modifiedModelAllLengths.xml");
		crf1.trainFeatures("test/input/multipleSequence/flyTwice");

		Conrad crf2 = new Conrad("test/input/multipleSequence/modifiedModelAllLengths2.xml");
		crf2.trainFeatures("test/input/multipleSequence/flyTwice");

		// Do an evaluation on a single sequence and save the alpha files
		crf1.trainWeights("test/input/multipleSequence/flyOnce");
		crf2.trainWeights("test/input/multipleSequence/flyTwice");

		// Do an evaluation on multiple sequences (where the first was the last) and save the alpha files
		String[] seq1 = FileUtil.readFile("test/working/multSeq1.txt").split("\n");
		String[] seq2 = FileUtil.readFile("test/working/multSeq2.txt").split("\n");

		int offset=0;
		while(seq2[offset].charAt(5)=='0') {
			++offset;
		}
		
		for(int i=0; i<seq1.length; ++i) {
			assertEquals(seq1[i].substring(6),seq2[i+offset].substring(6));
		}
	}
}
