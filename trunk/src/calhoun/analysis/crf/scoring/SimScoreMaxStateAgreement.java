package calhoun.analysis.crf.scoring;

import calhoun.analysis.crf.LocalPathSimilarityScore;
import calhoun.analysis.crf.io.TrainingSequence;

public class SimScoreMaxStateAgreement implements LocalPathSimilarityScore {

	public double evaluate(int yprev, int y, TrainingSequence<?> seq, int pos) {
		if (pos == 0) { return 0.0; }
		if (y == seq.getY(pos)) { return 1.0; }
		return 0.0;
	}

}
