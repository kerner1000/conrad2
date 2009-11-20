package calhoun.analysis.crf.test;

import java.util.List;

import junit.framework.TestCase;
import calhoun.analysis.crf.Conrad;
import calhoun.analysis.crf.features.interval13.StateTransitionsInterval13;
import calhoun.analysis.crf.io.TrainingSequence;

public class Interval13FeatureTest extends TestCase {

	public void testStateTransitions() throws Exception {
		Conrad cr = new Conrad("test/input/interval13/config/ssbaselineCPD.xml");
		List<? extends TrainingSequence<? extends Character>> data = (List<? extends TrainingSequence<? extends Character>>) cr.getInputHandler().readTrainingData("test/input/interval13/data/oneGeneTrain.interval13.txt");
		
		StateTransitionsInterval13 f = new StateTransitionsInterval13();
		f.train(0, cr.getModel(), data);
		assertEquals(Math.log(1/4.0), f.getEndProb());
		assertEquals(Math.log(3/4.0), f.getIntronProb());
	}
}
