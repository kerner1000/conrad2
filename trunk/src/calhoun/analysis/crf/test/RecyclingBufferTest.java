package calhoun.analysis.crf.test;

import calhoun.analysis.crf.solver.RecyclingBuffer;
import calhoun.util.AbstractTestCase;

public class RecyclingBufferTest extends AbstractTestCase {

	public void testRecyclingBuffer() throws Exception {
		RecyclingBuffer<Integer> i = new RecyclingBuffer<Integer>(new Integer[3]);
		i.addFirst(3);
		i.addFirst(2);
		assertEquals(2, (int) i.get(0));
		assertEquals(3, (int) i.get(1));
		i.addFirst(1);
		i.addFirst(0);
		assertEquals(0, (int) i.get(0));
		assertEquals(1, (int) i.get(1));
	}

	public void testRecyclingBufferWithArray() throws Exception {
		Integer[][] vals = new Integer[3][1];
		Integer[] val = new Integer[1];
		RecyclingBuffer<Integer[]> i = new RecyclingBuffer<Integer[]>(vals);
		val[0]=3;
		val = i.addFirst(val);
		val[0]=2;
		val = i.addFirst(val);
		assertEquals(2, i.get(0)[0].intValue());
		assertEquals(3, i.get(1)[0].intValue());
		val[0]=1;
		val = i.addFirst(val);
		val[0]=0;
		val = i.addFirst(val);
		assertEquals(0, i.get(0)[0].intValue());
		assertEquals(1, i.get(1)[0].intValue());
		val[0]=-1;
		val = i.addFirst(val);
		val[0]=-2;
		val = i.addFirst(val);
		assertEquals(-2, i.get(0)[0].intValue());
		assertEquals(-1, i.get(1)[0].intValue());
	}
}
