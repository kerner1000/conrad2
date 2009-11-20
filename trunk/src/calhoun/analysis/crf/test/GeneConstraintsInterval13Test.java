package calhoun.analysis.crf.test;

import java.util.ArrayList;
import java.util.List;

import calhoun.analysis.crf.Conrad;
import calhoun.analysis.crf.features.interval13.GeneConstraintsInterval13;
import calhoun.analysis.crf.io.StringInput;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.analysis.crf.solver.check.ArrayFeatureList;
import calhoun.util.AbstractTestCase;

/** Tests that CRF is working with valid probabilities - the sum of all possible labelings is 1.
 *  
 * Test that the code to walk through only the valid paths works correctly.
 * Uses a two state model that disallows transitions to self.  010101... or 101010... are the only allowed paths. */
public class GeneConstraintsInterval13Test extends AbstractTestCase {

	
	public void testGeneConstraintsTraining() throws Exception {

		Conrad crf = new Conrad("test/input/interval13/config/markov.xml");
		
		List<? extends TrainingSequence<?>> train1 =
			StringInput.prepareData(
				"000000002222222666661111100000000" + "\n" +
 				"ACACACACATGCACAGTCAGACACATAGACACA" + "\n" +
				"00000000077777CCCCCC7777000000000000" + "\n" +
				"ACACACTTACACACCTACACACATACACACACACAC" + "\n");
		
		System.out.println(train1);
		
		GeneConstraintsInterval13 gc = new GeneConstraintsInterval13();
		
		List<TrainingSequence<Character>> train2 = new ArrayList<TrainingSequence<Character>>();
		
		gc.train(0,crf.getModel(),train2);
	}

	
	public void testGeneConstraintsEvaluation() throws Exception {

		Conrad crf = new Conrad("test/input/interval13/config/markov.xml");
		
		GeneConstraintsInterval13 gc = new GeneConstraintsInterval13();
		List<? extends TrainingSequence<Character>> data = (List<? extends TrainingSequence<Character>>) crf.getInputHandler().readTrainingData("test/input/interval13/data/oneGeneTrain.interval13.txt");
		crf.trainFeatures(data);
		gc.train(0, crf.getModel(), data);
		
		ArrayFeatureList result = new ArrayFeatureList(crf.getModel());

		
		// Check that mod3 stuff done correctly for plus strand donor sites
		result.clear();
		gc.evaluateEdge(data.get(0).getInputSequence(), 3, 1, 4, result);
		assertTrue(result.isValid());

		result.clear();
		gc.evaluateEdge(data.get(0).getInputSequence(), 3, 2, 5, result);
		assertTrue(result.isValid());
		
		result.clear();
		gc.evaluateEdge(data.get(0).getInputSequence(), 3, 3, 6, result);
		assertTrue(result.isValid());
		
		result.clear();
		gc.evaluateEdge(data.get(0).getInputSequence(), 3, 1, 5, result);
		assertFalse(result.isValid());

		
		// check taht mod3 stuff done correctly for minus strand acceptor sites
		result.clear();
		gc.evaluateEdge(data.get(0).getInputSequence(), 18, 7, 10, result);
		assertTrue(result.isValid());
		
		result.clear();
		gc.evaluateEdge(data.get(0).getInputSequence(), 18, 8, 12, result);
		assertTrue(result.isValid());
		
		result.clear();
		gc.evaluateEdge(data.get(0).getInputSequence(), 18, 9, 11, result);
		assertTrue(result.isValid());
		
		result.clear();
		gc.evaluateEdge(data.get(0).getInputSequence(), 18, 7, 11, result);
		assertFalse(result.isValid());
		
		// check that some plus strand stop codons get invalidated for the exon state,
		// but only invalidated on third position and for exons of correct cut.
		result.clear();
		gc.evaluateNode(data.get(0).getInputSequence(), 7, 3, result);
		assertFalse(result.isValid());		

		result.clear();
		gc.evaluateNode(data.get(0).getInputSequence(), 7, 2, result);
		assertTrue(result.isValid());		
		
		result.clear();
		gc.evaluateNode(data.get(0).getInputSequence(), 6, 3, result);
		assertTrue(result.isValid());		

		result.clear();
		gc.evaluateEdge(data.get(0).getInputSequence(), 7, 3, 3, result);
		assertFalse(result.isValid());	

		result.clear();
		gc.evaluateEdge(data.get(0).getInputSequence(), 7, 4, 4, result);
		assertTrue(result.isValid());	
		
		// check that some minus strand stop codons get invalidated for the exon state,
		// but only invalidated on third position and for exons of correct cut.
		result.clear();
		gc.evaluateNode(data.get(0).getInputSequence(), 4, 8, result);
		assertFalse(result.isValid());		

		result.clear();
		gc.evaluateNode(data.get(0).getInputSequence(), 4, 7, result);
		assertTrue(result.isValid());		
		
	}
}
