package calhoun.analysis.crf.test;

import java.util.List;

import junit.framework.TestCase;
import calhoun.analysis.crf.Conrad;
import calhoun.analysis.crf.features.interval29.StateTransitionsInterval29;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.analysis.crf.io.SequenceConverter;

public class Interval29FeatureTest extends TestCase {

	public void testStateTransitions() throws Exception {
		Conrad cr = new Conrad("test/input/interval29/config/ssbaselineCPD2Int29.xml");
		List<? extends TrainingSequence<? extends Character>> data = (List<? extends TrainingSequence<? extends Character>>) cr.getInputHandler().readTrainingData("test/input/interval13/data/oneGeneTrain.interval13.txt");
		for (TrainingSequence<? extends Character> seq : data) {
			seq.setY(SequenceConverter.convertSeqFromInterval13ToInterval29(seq.getY()));
		}
		StateTransitionsInterval29 f = new StateTransitionsInterval29();
		f.train(0, cr.getModel(), data);
		assertEquals(Math.log(1/4.0), f.getEndProb());
		assertEquals(Math.log(3/4.0), f.getIntronProb());
	}
}
