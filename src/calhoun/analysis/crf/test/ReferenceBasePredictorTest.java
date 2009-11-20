package calhoun.analysis.crf.test;

import java.util.List;

import calhoun.analysis.crf.Conrad;
import calhoun.analysis.crf.features.interval13.ReferenceBasePredictorInterval13;
import calhoun.analysis.crf.features.supporting.LogProbLookup;
import calhoun.analysis.crf.io.InputSequenceCharacter;
import calhoun.analysis.crf.io.StringInput;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.util.AbstractTestCase;

/** Tests that CRF is working with valid probabilities - the sum of all possible labelings is 1.
 *  
 * Test that the code to walk through only the valid paths works correctly.
 * Uses a two state model that disallows transitions to self.  010101... or 101010... are the only allowed paths. */
public class ReferenceBasePredictorTest extends AbstractTestCase {


	public void testLogProbLookupEasy() throws Exception {
		
		InputSequenceCharacter seq = new InputSequenceCharacter("CTCTCTCTCTCTCTCTCTC");

		// lp1 will measure single nucleotide frequencies for the positive strand
		// in this sequence there are 10 C's and 9 T's, but with pseudocounts of
		// 1.0 this bumps up to 11 C's 10 T's 1 A and 1 G, total 23 letters.
		LogProbLookup lp1 = new LogProbLookup(0,1.0);	
		for (int pos=0; pos<seq.length(); pos++) {
			lp1.increment(seq,pos,true);
		}
		lp1.finalize();

		assertEquals(lp1.lookup(seq,0,true),-0.737598,0.001);  // log of 11/23 = -0.737599 
		assertEquals(lp1.lookup(seq,0,false),-3.135494,0.001); // log of 1/23 is -3.1354
		assertEquals(lp1.lookup(seq,1,true),-0.832909,0.001);  // log of 10/23 is -0.832909
		assertEquals(lp1.lookup(seq,1,false),-3.135494,0.001); // log of 1/23 is -3.1354
		

		// lp2 is like lp1 except the pseudocounts will be 0.25 instead of 1.0
		LogProbLookup lp2 = new LogProbLookup(0,0.25);	
		for (int pos=0; pos<seq.length(); pos++) {
			lp2.increment(seq,pos,true);
		}
		lp2.finalize();

		assertEquals(lp2.lookup(seq,0,true),-0.66845,0.001);  // log of C 10.25/20
		assertEquals(lp2.lookup(seq,0,false),-4.38203,0.001); // log of G 0.25/20
		assertEquals(lp2.lookup(seq,1,true),-0.77111,0.001);  // log of T 9.25/20
		assertEquals(lp2.lookup(seq,1,false),-4.38203,0.001); // log of A 0.25/20 	
	}
	
	
	public void testLogProbLookupHard() throws Exception {
		
		InputSequenceCharacter seq1 = new InputSequenceCharacter("ACGTNCGTGTTCCATGGTAAC");
		InputSequenceCharacter seq2 = new InputSequenceCharacter("GNTTACA");
			
		
		// lp1 will measure probability of a letter based on previous two, pseudocounts=1.0
		// lp1 being trained on both strands of seq1
		LogProbLookup lp1 = new LogProbLookup(2,1.0);	
		for (int pos=0; pos<seq1.length(); pos++) {
			lp1.increment(seq1,pos,true);
			lp1.increment(seq1,pos,false);
		}
		lp1.finalize();

		assertEquals(lp1.lookup(seq2,0,false),0.0,0.001); // should pop up as missing data, default to 0.0
		assertEquals(lp1.lookup(seq2,3,true),0.0,0.001); // should pop up as missing data b/c of N at position 1, default to 0.0
		
		assertEquals(lp1.lookup(seq2,4,true),-1.09861,0.001);
		//the history here is "TT".  In the training sequence occurs once on the positive strand,
		// yielding C, and once on the negative strand yeiding the revcomp(A).  
		// Thus the pseudocounts are: C:2, A:2, G:1, T:1
		// Here it's A so I expect log(2/6) = -1.09861
	}	
	
	
	public void testReferenceBasePredictor() throws Exception {
	
		String configFile = "test/input/interval13/config/markov.xml";
		
		Conrad crf = new Conrad(configFile);
		
		List<? extends TrainingSequence<Character>> train1 =
			StringInput.prepareData(
				"000000002222222666661111100000000" + "\n" +
 				"ACACACACATGCACAGTCAGACACATAGACACA" + "\n" +
				"00000000077777CCCCCC7777000000000000" + "\n" +
				"ACACACTTACACACCTACACACATACACACACACAC" + "\n");
		
		System.out.println(train1);
		
		ReferenceBasePredictorInterval13 bp = new ReferenceBasePredictorInterval13();
		
		bp.train(0,crf.getModel(),train1);
		
	}
	
}
