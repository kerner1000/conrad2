package calhoun.analysis.crf.test.slow;

import java.util.ArrayList;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.features.supporting.MaxentMotifModel;
import calhoun.analysis.crf.test.CRFIOTest;
import calhoun.util.AbstractTestCase;

public class SlowTestMaxentMotifModel extends AbstractTestCase {
	private static final Log log = LogFactory.getLog(CRFIOTest.class);
	boolean debug = log.isDebugEnabled();
	
	public void testSlowMaxentModel() throws Exception {
		ArrayList<int[]> motifExamples = new ArrayList<int[]>();
		motifExamples.add(new int[]{0,0,0,0,0,0,0,0,0,0,0});
		//motifExamples.add(new int[]{1,1,0});
		//motifExamples.add(new int[]{1,0,1});
		//motifExamples.add(new int[]{0,1,1});
		motifExamples.add(new int[]{1,1,1,1,1,1,1,1,1,1,1});
		
		double[] pp = MaxentMotifModel.trainMaxentDistributionUsingAllPairwiseConstraints(motifExamples,11,30,1.0);
		
		assertEquals(pp.length,4194304);	
	}

}
