package calhoun.analysis.crf.features.tricycle13;

import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.AbstractFeatureManager;
import calhoun.analysis.crf.CacheStrategySpec;
import calhoun.analysis.crf.FeatureList;
import calhoun.analysis.crf.FeatureManagerEdge;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.CacheStrategySpec.CacheStrategy;
import calhoun.analysis.crf.features.supporting.MarkovPredictorLogprob;
import calhoun.analysis.crf.features.tricycle13.EmissionMarkovFeature.MarkovHistory;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.seq.KmerHasher;
import calhoun.util.Assert;

public class CodingStopFeature extends AbstractFeatureManager<Character> implements FeatureManagerEdge<Character> {
	private static final long serialVersionUID = -7659288739348604129L;
	private static final Log log = LogFactory.getLog(CodingStopFeature.class);
	boolean debug = log.isDebugEnabled();
	
	int startIx;  // The index of the first feature managed by this FeatureManager
	ModelManager model;
	
	KmerHasher   h; // for a single letter
	
	float[][] pluslogprob; /* logprob[i][j] is log probability of base j at position i, i=0..(span-1), j=0..3. */
	float[][] minuslogprob;

	int[][] pc; // for plus corrections.  Second dimension is of length 2; First entry is where to look relative to the transition; second entry is what the hidden state in that position should be.  
	int[][] mc; // minus corrections

	int stateIntergenic;
	int stateExon3;
	int stateExon3m;

	MarkovPredictorLogprob predictorlp;
	
	///////////////////////////////////// Class variables above, methods below //////////

	
	public CodingStopFeature() {
	}

	public void setHistory(MarkovHistory markovHistory) {
		this.predictorlp = new MarkovPredictorLogprob(markovHistory.convert());
	}
	
	public CodingStopFeature( List<int[]> markovhistory ) {
		this.predictorlp = new MarkovPredictorLogprob(markovhistory);	
	}


	public int getNumFeatures() {
		return 1;
	}	
	
	
	public String getFeatureName(int featureIndex) {
		Assert.a(featureIndex == startIx);
		return "CodingStopFeature";
	}

	
	
	public void evaluateEdge(InputSequence<? extends Character> seq, int ix, int prevState, int state, FeatureList result) {
		int len = seq.length();
		if(prevState == stateExon3 && state == stateIntergenic) {
			// first the positive strand
			float pval = 0;
			if ((ix>=3) && (ix<len-5) ) {

				char a = seq.getX(ix-3);
				char b = seq.getX(ix-2);
				char c = seq.getX(ix-1);
				if ( (a=='T') && (b=='A') && (c=='G') ) { 	pval = (float) Math.log(0.33333333); }
				else if ( (a=='T') && (b=='G') && (c=='A') ) { pval = (float) Math.log(0.33333333); }
				else if ( (a=='T') && (b=='A') && (c=='A') ) { pval = (float) Math.log(0.33333333); }
				else { pval = (float) -50.0; }
			
				for (int i=0; i<6; i++) {
					char d = seq.getX(ix+i);
					pval += pluslogprob[i][h.hash(d)];
				}
				
				for (int i=0; i<9; i++) {
					pval -= predictorlp.logprob(pc[i][1],seq,ix-pc[i][0]);					
				}
			}
			result.addFeature(startIx, pval);
		}
		else if(state == stateExon3m && prevState == stateIntergenic) {
			// now the negative strand;
			float mval = 0;
			if ((ix>=6) && (ix<len-2) ) {

				char a = seq.getX(ix+0);
				char b = seq.getX(ix+1);
				char c = seq.getX(ix+2);
				if ( (a=='T') && (b=='T') && (c=='A') ) { 	mval = (float) Math.log(0.33333333); }
				else if ( (a=='C') && (b=='T') && (c=='A') ) { mval = (float) Math.log(0.33333333); }
				else if ( (a=='T') && (b=='C') && (c=='A') ) { mval = (float) Math.log(0.33333333); }
				else { mval = (float) -50.0; }
			
				for (int i=0; i<6; i++) {
					char d = seq.getX(ix+i-6);
					mval += minuslogprob[i][h.hash(d)];
				}
				
				for (int i=0; i<9; i++) {
					mval -= predictorlp.logprob(mc[i][1],seq,ix-mc[i][0]);					
				}
			}
			result.addFeature(startIx, mval);
		}
	}

	
	public void train(int startingIndex, ModelManager modelInfo, List<? extends TrainingSequence<? extends Character>> data) {

		startIx = startingIndex;
		model = modelInfo;

		stateIntergenic = model.getStateIndex("intergenic");
		stateExon3 =  model.getStateIndex("exon3");
		stateExon3m =  model.getStateIndex("exon3m");

		predictorlp.train(data);
		h = new KmerHasher(KmerHasher.ACGTN, 1);		
		
		
		/*
		# stop signal positive strand TAG|xxxxxx   
		L 9 3 exon3 intergenic
		# stop signal negative strand xxxxxx|CTA
		L 9 6 intergenic exon3m

		# stop signal positive strand TAG|xxxxxx  
		L exon1 exon2 exon3 intergenic intergenic intergenic intergenic intergenic intergenic
		# stop signal negative strand xxxxxx|CTA
		L intergenic intergenic intergenic intergenic intergenic intergenic exon3m exon2m exon1m
		*/
		
		pc = new int[9][2];
		pc[0][0] = -3;  // pc[0][1] = model.getStateIndex("exon1");
		pc[1][0] = -2;   pc[1][1] = model.getStateIndex("exon2");
		pc[2][0] = -1;   pc[2][1] = model.getStateIndex("exon3");
		pc[3][0] =  0;   pc[3][1] = model.getStateIndex("intergenic");
		pc[4][0] =  1;   pc[4][1] = model.getStateIndex("intergenic");
		pc[5][0] =  2;   pc[5][1] = model.getStateIndex("intergenic");
		pc[6][0] =  3;   pc[6][1] = model.getStateIndex("intergenic");
		pc[7][0] =  4;   pc[7][1] = model.getStateIndex("intergenic");
		pc[8][0] =  5;   pc[8][1] = model.getStateIndex("intergenic");

		mc = new int[9][2];
		mc[0][0] = -6;   mc[0][1] = model.getStateIndex("intergenic");
		mc[1][0] = -5;   mc[1][1] = model.getStateIndex("intergenic");
		mc[2][0] = -4;   mc[2][1] = model.getStateIndex("intergenic");
		mc[3][0] = -3;   mc[3][1] = model.getStateIndex("intergenic");
		mc[4][0] = -2;   mc[4][1] = model.getStateIndex("intergenic");
		mc[5][0] = -1;   mc[5][1] = model.getStateIndex("intergenic");
		mc[6][0] =  0;   mc[6][1] = model.getStateIndex("exon3m");
		mc[7][0] =  1;   mc[7][1] = model.getStateIndex("exon2m");
		mc[8][0] =  2;   mc[8][1] = model.getStateIndex("exon1m");
			
		pluslogprob = new float[6][h.range()];
		minuslogprob = new float[6][h.range()];		
		
		for (int i=0; i<6; i++) {
			for (int j=0; j<h.range(); j++) {
				pluslogprob[i][j] = (float) 1.0;
				minuslogprob[i][j] = (float) 1.0;				
			}
		}
		
		

		// In English, what I want to do is this.  Loop through all of the training data, once for each Feature.
		// While so doing, look for any positions where one of the allowed transitions for that feature occurs.
		// At such positions, increment the counts for logprob.
		for(TrainingSequence<? extends Character> seq : data) {
			int len = seq.length();		
			
			for (int ix=3; ix<(len-5); ix++) {
				int yprev = seq.getY(ix-1);
				int y = seq.getY(ix);
				if (  (yprev == model.getStateIndex("exon3") ) && (y == model.getStateIndex("intergenic") )  ) {	
					for (int k=0; k<6; k++) {
						char c = seq.getX(ix + k);
						pluslogprob[k][h.hash(c)] += 1.0;
					}
				}
			}
			
			for (int ix=6; ix<(len-3); ix++) {
				int yprev = seq.getY(ix-1);
				int y = seq.getY(ix);
				if (  (yprev == model.getStateIndex("intergenic") ) && (y == model.getStateIndex("exon3m") )  ) {	
					for (int k=0; k<6; k++) {
						char c = seq.getX(ix + k - 6);
						minuslogprob[k][h.hash(c)] += 1.0;
					}
				}
			}	
				
		}

		for (int k=0; k<6; k++) {
			float totalp = 0;
			for (int j=0; j<h.range(); j++) { totalp += pluslogprob[k][j]; }
			Assert.a(totalp > 0);
			for (int j=0; j<h.range(); j++) { 
				pluslogprob[k][j] = (float) (Math.log(pluslogprob[k][j]) - Math.log(totalp));
			}
		}
		
		for (int k=0; k<6; k++) {			
			float totalm = 0;
			for (int j=0; j<h.range(); j++) { totalm += minuslogprob[k][j]; }
			Assert.a(totalm > 0);
			for (int j=0; j<h.range(); j++) { 
				minuslogprob[k][j] = (float) (Math.log(minuslogprob[k][j]) - Math.log(totalm));
			}	
		}		
	}
	@Override
	public CacheStrategySpec getCacheStrategy() {
		return new CacheStrategySpec(CacheStrategy.UNSPECIFIED);
	}
	
}

