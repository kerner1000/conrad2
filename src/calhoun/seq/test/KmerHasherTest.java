package calhoun.seq.test;

import calhoun.seq.KmerHasher;
import calhoun.util.AbstractTestCase;

public class KmerHasherTest extends AbstractTestCase  {
	public void testDnaLengthOne() {
		KmerHasher h = new KmerHasher(KmerHasher.DNA, 1);
		assertEquals(0, h.hash("a", 0));
		assertEquals(0, h.hash("a", 0));
		assertEquals(1, h.hash("C", 0));
		assertEquals(2, h.hash("G", 0));
		assertEquals(3, h.hash("T", 0));
	}

	public void testLetterLengthOne() {
		KmerHasher h = new KmerHasher(KmerHasher.LETTERS, 1);
		assertEquals(0, h.hash("a", 0));
		assertEquals(0, h.hash("A", 0));
		assertEquals(1, h.hash("b", 0));
		assertEquals(25, h.hash("Z", 0));
	}

	public void testDnaLengthTwo() {
		KmerHasher h = new KmerHasher(KmerHasher.DNA, 2);
		assertEquals(0, h.hash("aa", 0));
		assertEquals(1, h.hash("ac", 0));
		assertEquals(4, h.hash("ca", 0));
		assertEquals(5, h.hash("cc", 0));
		assertEquals(15, h.hash("TT", 0));
	}

	public void testLettersLengthTwo() {
		KmerHasher h = new KmerHasher(KmerHasher.LETTERS, 2);
		assertEquals(0, h.hash("aa", 0));
		assertEquals(1, h.hash("ab", 0));
		assertEquals(26, h.hash("ba", 0));
		assertEquals(26*25+1, h.hash("zb", 0));
	}

	public void testCallTypes() {
		KmerHasher h = new KmerHasher(KmerHasher.DNA, 4);
		assertEquals(130, h.hash("gaag", 0));
		assertEquals(130, h.hash("aagaagaa", 2));
		assertEquals(130, h.hash(new char[] {'g', 'a', 'a', 'g'}));
	}

	public void testShift() {
		KmerHasher h = new KmerHasher(KmerHasher.DNA, 4);
		// gaag hashes to 130
		int hash = h.hash("agaa", 0);
		assertEquals(32, hash);
		assertEquals(130, h.shiftHash('g', hash));
		assertEquals(hash, h.reverseShiftHash('a', 130));
	}

	public void testReverse() {
		/*KmerHasher h = new KmerHasher(KmerHasher.LETTERS, 4);
		// gaag hashes to 130
		int hash = h.hash("bb", 0);
		assertEquals(27, hash);
		assertEquals(130, h.shiftHash('g', hash));
		assertEquals(hash, h.reverseShiftHash('a', 130));*/
	}
}
