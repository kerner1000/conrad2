package calhoun.analysis.crf.scoring;

import calhoun.analysis.crf.LocalPathSimilarityScore;
import calhoun.analysis.crf.io.TrainingSequence;

public class SimScoreMinCodingMiscallsSV13 implements LocalPathSimilarityScore {

	public double evaluate(int yprev, int y, TrainingSequence<?> seq, int pos) {
		if (pos == 0) { return 0.0; }
		int realy = seq.getY(pos);
		if (isCodingPlus(y) && (!isCodingPlus(realy))) { return -1.0; }
		if ((!isCodingPlus(y)) && isCodingPlus(realy)) { return -1.0; }
		if (isCodingMinus(y) && (!isCodingMinus(realy))) { return -1.0; }
		if ((!isCodingMinus(y)) && isCodingPlus(realy)) { return -1.0; }
		return 0.0;
	}

	private boolean isCodingMinus(int y) {
		if (y==7) { return true; }
		if (y==8) { return true; }
		if (y==9) { return true; }
		return false;
	}

	private boolean isCodingPlus(int y) {
		if (y==1) { return true; }
		if (y==2) { return true; }
		if (y==3) { return true; }
		return false;
	}

}
	
