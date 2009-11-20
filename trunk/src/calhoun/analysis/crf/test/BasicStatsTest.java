package calhoun.analysis.crf.test;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.statistics.BasicStats;
import calhoun.util.AbstractTestCase;

public class BasicStatsTest extends AbstractTestCase {
	private static final Log log = LogFactory.getLog(CRFIOTest.class);
	boolean debug = log.isDebugEnabled();
	
	public void testMean() throws Exception {
		double[] x = new double[]{1,2,3,4,5};
		assertEquals(BasicStats.meanDoubleArray(x),3,0.000001);
		assertEquals(Math.pow(4,3),64,0.000001);
	}
	
	public void testMedian() throws Exception {
		double[] x = new double[]{1,2,3,4,50};
		assertEquals(BasicStats.medianDoubleArray(x),3,0.000001);
	}

}
