package calhoun.analysis.crf.test;

import calhoun.analysis.crf.statistics.PredictedActualBinaryContingencyTable;
import calhoun.util.AbstractTestCase;

public class PredictedActualBinaryContingencyTableTest extends AbstractTestCase {

	
	public void testFull() throws Exception {
		PredictedActualBinaryContingencyTable box = new PredictedActualBinaryContingencyTable();
		box.set(100,5,10,1000);
		System.out.println(box.summarize());	
		
		PredictedActualBinaryContingencyTable box2 = new PredictedActualBinaryContingencyTable();
		box2.set(1000,5,10,100);
		System.out.println(box2.summarize());	
	}
	
	public void testPartial() throws Exception {
		PredictedActualBinaryContingencyTable box = new PredictedActualBinaryContingencyTable();
		box.set(100,5,10);
		System.out.println(box.summarize());	
		
		PredictedActualBinaryContingencyTable box2 = new PredictedActualBinaryContingencyTable();
		box2.set(1000,5,10);
		System.out.println(box2.summarize());	
	}	

	public void testIncremental() throws Exception {
		PredictedActualBinaryContingencyTable box = new PredictedActualBinaryContingencyTable();
		for (int i=0; i<100; i++) { box.incrementTP(); }
		for (int i=0; i<5; i++) { box.incrementFP(); }
		for (int i=0; i<10; i++) { box.incrementFN(); }
		for (int i=0; i<1000; i++) { box.incrementTN(); }
		box.freeze();
		System.out.println(box.summarize());		
	}
	
	public void testIncrementalPartial() throws Exception {
		PredictedActualBinaryContingencyTable box = new PredictedActualBinaryContingencyTable();
		box.forgetTN();
		for (int i=0; i<100; i++) { box.incrementTP(); }
		for (int i=0; i<5; i++) { box.incrementFP(); }
		for (int i=0; i<10; i++) { box.incrementFN(); }
		box.freeze();
		System.out.println(box.summarize());				
	}

}
