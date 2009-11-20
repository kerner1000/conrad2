package calhoun.analysis.crf.test;

import java.util.ArrayList;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.features.supporting.MaxentMotifModel;
import calhoun.analysis.crf.statistics.BasicStats;
import calhoun.util.AbstractTestCase;
import calhoun.util.FileUtil;

public class MaxentMotifModelTest extends AbstractTestCase {
	private static final Log log = LogFactory.getLog(CRFIOTest.class);
	boolean debug = log.isDebugEnabled();
	
	public void testMaxentModel() throws Exception {
		ArrayList<int[]> motifExamples = new ArrayList<int[]>();
		motifExamples.add(new int[]{0,0,0});
		//motifExamples.add(new int[]{1,1,0});
		//motifExamples.add(new int[]{1,0,1});
		//motifExamples.add(new int[]{0,1,1});
		motifExamples.add(new int[]{1,1,1});
		
		double[] pp = MaxentMotifModel.trainMaxentDistributionUsingAllPairwiseConstraints(motifExamples,3,30,1.0);
		
		assertEquals(pp.length,64);
		
		for (int a=0; a<4; a++) {
			for (int b=0; b<4; b++) {
				for (int c=0; c<4; c++) {
					int h = 16*a+4*b+c;
					if (pp[h] > 0.001) {
						System.out.println("" + a + " " + b + " " + c + "   " + pp[h]);
					}
				}
			}
		}
		
		ArrayList<int[]> motifExamples2 = new ArrayList<int[]>();
		
		double[] qq = MaxentMotifModel.trainMaxentDistributionUsingAllPairwiseConstraints(motifExamples2,3,30,1.0);
		
		assertEquals(qq.length,64);
		
		for (int a=0; a<4; a++) {
			for (int b=0; b<4; b++) {
				for (int c=0; c<4; c++) {
					int h = 16*a+4*b+c;
					if (qq[h] > 0.001) {
						System.out.println("" + a + " " + b + " " + c + "   " + qq[h]);
					}
				}
			}	
		}
	}
	
	public void testMaxentDonorModel() throws Exception {
		// We begin with lots of examples of donor sites, build the maximum entropy distribution
		// based on all pairwise marginals using Java, and compare results with the same
		// maxent distribution (of length 16384 = 4^7) using Matlab.

		String fileName ="test/input/donor_examples.txt";
		String[][] preDonorExamples = FileUtil.readFlatFile(fileName);
		int nExamples = preDonorExamples.length;
		assertEquals(nExamples,1098);
		int span = 7;
		System.out.println("Number of donor examples is " + nExamples);

		ArrayList<int[]> donorExamples = new ArrayList<int[]>();		
		for (int j=0; j<nExamples; j++) {
			assertEquals(preDonorExamples[j].length,span);
			int[] example = new int[span];
			for (int k=0; k<span; k++) {
				example[k] = (int) Math.round(Double.parseDouble(preDonorExamples[j][k])) - 1;
			}
			donorExamples.add(example);
		}
		
		
		
		String fileName2 = "test/input/donor_maxent.txt";
		double[] x = FileUtil.readDoublesFromSingleTabbedLine(fileName2); 
		System.out.println("Length is " + x.length  +"     Sum is " + BasicStats.sumDoubleArray(x));
		
		double[] y = MaxentMotifModel.trainMaxentDistributionUsingAllPairwiseConstraints(donorExamples,span,200,0.0);
		
		double diff = BasicStats.L1Distance(x,y);
		System.out.println("The L1 difference between the Matlab and Java maxent distributions was " + diff);
		assert(diff < 0.01);
		
		System.out.println("The max and argmax of the Java distribution are " + BasicStats.max(y) + "    " + BasicStats.argmax(y));
		System.out.println("The max and argmax of the Matlab distribution are " + BasicStats.max(x) + "    " + BasicStats.argmax(x));
		
	}
	
	

}
