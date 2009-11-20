package calhoun.analysis.crf.test;

import calhoun.analysis.crf.features.supporting.LogProbLookup;
import calhoun.analysis.crf.io.InputSequenceCharacter;
import calhoun.seq.KmerHasher;
import calhoun.util.AbstractTestCase;

/** Tests that CRF is working with valid probabilities - the sum of all possible labelings is 1.
 *  
 * Test that the code to walk through only the valid paths works correctly.
 * Uses a two state model that disallows transitions to self.  010101... or 101010... are the only allowed paths. */
public class LogprobLookupTest extends AbstractTestCase {

	
	public void testLogProbLookupRC() throws Exception {
		String forward = "TGTTGGTACGCTTGCGGCTCTGCTGCAGCGAAAAAAAAGATCGAAATGACCAG";
		String reverse = "CTGGTCATTTCGATCTTTTTTTTCGCTGCAGCAGAGCCGCAAGCGTACCAACA";
		String computedReverse = KmerHasher.reverseComplement(forward);

		int len = forward.length();		
		assertEquals(reverse,computedReverse);
		

		
		InputSequenceCharacter forwardSeq = new InputSequenceCharacter(forward);
		InputSequenceCharacter reverseSeq = new InputSequenceCharacter(reverse);		
	
		LogProbLookup lf = new LogProbLookup(2,1.0);

		LogProbLookup lr = new LogProbLookup(2,1.0);
		
		for (int j=0; j<len; j++) {
			lf.increment(forwardSeq,j,true);
			lr.increment(reverseSeq,j,false);
		}
		
		lf.finalize(); 
		lr.finalize();
		
		for (int j=0; j<len; j++) {
			assertEquals(lf.lookup(forwardSeq,j,true),lr.lookup(forwardSeq,j,true),0.0001);
			assertEquals(lf.lookup(forwardSeq,j,false),lr.lookup(forwardSeq,j,false),0.0001);			
			assertEquals(lf.lookup(forwardSeq,j,true),lr.lookup(reverseSeq,len-j-1,false),0.0001);
		}
		
		
		System.out.println("hello");
		
	}	
	
}
