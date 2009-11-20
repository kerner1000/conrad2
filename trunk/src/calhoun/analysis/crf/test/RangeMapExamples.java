package calhoun.analysis.crf.test;

import calhoun.analysis.crf.io.IntervalInputSequence.IntervalRangeMapValue;
import calhoun.util.AbstractTestCase;
import calhoun.util.RangeMap;

public class RangeMapExamples extends AbstractTestCase {
	
	public void testRangeMapPositionQueryExample() throws Exception {
		RangeMap rm = new RangeMap();
		
		int[] starts = new int[]{5,10,15,20,12};
		int[] ends = new int[]{7,9,17,22,23};
		assertEquals(starts.length,ends.length);
		int n = starts.length;
		
		
		for (int j=0; j<n; j++) {
			IntervalRangeMapValue irmv = new IntervalRangeMapValue(5,10,1.0);
			rm.add(irmv.start,irmv.end,irmv);
		}
		
		
		
		
		
	}

}
