package calhoun.analysis.crf.scoring;

import calhoun.analysis.crf.LocalPathSimilarityScore;
import calhoun.analysis.crf.io.TrainingSequence;

public class SimScoreStateAndExonBoundariesInt13 implements LocalPathSimilarityScore {
	
	double stateBonus = 1.0;
	double boundaryPenalty = 50;
	
	public double evaluate(int yprev, int y, TrainingSequence<?> seq, int pos) {
		if (pos == 0) { return 0.0; }
		int realy = seq.getY(pos);
		int realyprev = seq.getY(pos-1);
		double score = y == realy ? stateBonus : 0.0; 
		
		if ((isPStart(yprev,y) ^ isPStart(realyprev,realy)) ||
				(isPDon(yprev,y)   ^ isPDon(realyprev,realy))	||
				(isPAcc(yprev,y)   ^ isPAcc(realyprev,realy))	||
				(isPStop(yprev,y)  ^ isPStop(realyprev,realy))	||
				(isMStart(yprev,y) ^ isMStart(realyprev,realy))	||
				(isMDon(yprev,y)   ^ isMDon(realyprev,realy))	||
				(isMAcc(yprev,y)   ^ isMAcc(realyprev,realy))	||
				(isMStop(yprev,y)  ^ isMStop(realyprev,realy))
				) { 
			score -= boundaryPenalty; 
		}
	
		return score;
	}

	private boolean isPStart(int yprev, int y) {
		if ( (yprev==0) && (y==1)) { return true; }
		return false;
	}

	private boolean isPDon(int yprev, int y) {
		if ( (yprev==1) && (y==4)) { return true; }
		if ( (yprev==2) && (y==5)) { return true; }
		if ( (yprev==3) && (y==6)) { return true; }
		return false;
	}

	private boolean isPAcc(int yprev, int y) {
		if ( (yprev==4) && (y==2)) { return true; }
		if ( (yprev==5) && (y==3)) { return true; }
		if ( (yprev==6) && (y==1)) { return true; }
		return false;
	}

	private boolean isPStop(int yprev, int y) {
		if ( (yprev==3) && (y==0)) { return true; }
		return false;
	}
	
	private boolean isMStart(int yprev, int y) {
		if ( (yprev==7) && (y==0)) { return true; }
		return false;
	}

	private boolean isMDon(int yprev, int y) {
		if ( (yprev==11) && (y==8)) { return true; }
		if ( (yprev==10) && (y==7)) { return true; }
		if ( (yprev==12) && (y==9)) { return true; }
		return false;
	}
	
	private boolean isMAcc(int yprev, int y) {
		if ( (yprev==9) && (y==11)) { return true; }
		if ( (yprev==8) && (y==10)) { return true; }
		if ( (yprev==7) && (y==12)) { return true; }
		return false;
	}
	
	private boolean isMStop(int yprev, int y) {
		if ( (yprev==0) && (y==9)) { return true; }
		return false;
	}
}
