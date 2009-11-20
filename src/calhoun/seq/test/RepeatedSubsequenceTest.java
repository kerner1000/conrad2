package calhoun.seq.test;

import java.util.Map;

import calhoun.seq.RepeatedSubsequence;
import calhoun.util.AbstractTestCase;

public class RepeatedSubsequenceTest extends AbstractTestCase  {
	public void testRepeatedSubsequence() {
		assertEquals("ATGC", RepeatedSubsequence.calc("ATGC"));
		assertEquals("A", RepeatedSubsequence.calc("AA"));
		assertEquals("A", RepeatedSubsequence.calc("AAA"));
		assertEquals("A", RepeatedSubsequence.calc("AAAA"));
		assertEquals("A", RepeatedSubsequence.calc("AAAAA"));
		assertEquals("A", RepeatedSubsequence.calc("AAAAAA"));
		assertEquals("AT", RepeatedSubsequence.calc("ATATATATAT"));
		assertEquals("ATATCG", RepeatedSubsequence.calc("ATATCGATATCG"));
		assertEquals("AT", RepeatedSubsequence.calc("ATATATATATAT"));
	}
	
	public void testRollingMatch() {
		assertEquals(true, RepeatedSubsequence.isRollingMatch("ATGC", "ATGC"));
		assertEquals(true, RepeatedSubsequence.isRollingMatch("CATG", "ATGC"));
		assertEquals(true, RepeatedSubsequence.isRollingMatch("GCAT", "ATGC"));
		assertEquals(true, RepeatedSubsequence.isRollingMatch("TGCA", "ATGC"));
		assertEquals(false, RepeatedSubsequence.isRollingMatch("TCGA", "ATGC"));
	}

	private String getCanonical(Map<String, Integer> m, String seq) {
		// Get the smallest repeated subsequence from this sequence
		seq = RepeatedSubsequence.calc(seq);

		// Check this subseqeunce vs. all others to see if we have a rolling match.
		for(String key : m.keySet()) {
			if(RepeatedSubsequence.isRollingMatch(key, seq)) {
				m.put(key, m.get(key) + 1);
				return key;
			}
		}

		// If no match, add it into the map
		m.put(seq, 1);
		return seq;
	}
}
