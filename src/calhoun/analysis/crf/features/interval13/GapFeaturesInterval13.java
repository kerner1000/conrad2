package calhoun.analysis.crf.features.interval13;

import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.AbstractFeatureManager;
import calhoun.analysis.crf.CacheStrategySpec;
import calhoun.analysis.crf.FeatureList;
import calhoun.analysis.crf.FeatureManagerNode;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.CacheStrategySpec.CacheStrategy;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.analysis.crf.io.MultipleAlignmentInputSequence;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.analysis.crf.io.MultipleAlignmentInputSequence.MultipleAlignmentColumn;
import calhoun.seq.KmerHasher;
import calhoun.util.Assert;

public class GapFeaturesInterval13 extends AbstractFeatureManager<MultipleAlignmentColumn> implements FeatureManagerNode<MultipleAlignmentColumn> {
	private static final long serialVersionUID = -7659288739348604129L;
	private static final Log log = LogFactory.getLog(GapFeaturesInterval13.class);
	boolean debug = log.isDebugEnabled();
	String[] featureNames = new String[] { "Frameshift coding", "Frameshift intron", "Frameshift exon", "Mod 3 gap coding", "Mod 3 gap intron", "Mod 3 gap exon"};
	
	/* Let G(i) and F(i) be booleans, one per position, that indicate whether there is
	 * a gap in the multiple alignment for which either 1) The first non-gap character
	 * of the reference sequence to the right of the gap is at position i, or 2) The last
	 * non-gap character of the reference sequence to the left of the gap is at position i;
	 * and 3) the gap is a multiple of 3.  G stands for "gap"
	 * 
	 * F(i) is the same thing except for a non-multiple of three length.
	 * F stands for "frameshifter".
	 * 
	 * The we define indicator features for the following conjunctions:
	 * 1) F(i) & (coding)
	 * 2) F(i) & (intronic)
	 * 3) F(i) & (intergenic)
	 * 4) G(i) & (coding)
	 * 5) G(i) & (intronic)
	 * 6) G(i) & (intergenic)
	 */
	
	// This code is essentially unchanged from the tricycle13 feature of similar name.
	
	int startIx;  
	ModelManager model;
	KmerHasher h = new KmerHasher(KmerHasher.ACGTN,1);
	
	int maxSeqLength;
	
	int nFeatures = 6;
	int nStates;
	
	Boolean[] gapboundary, frameshifter;
	Boolean[] isStateCoding, isStateIntronic, isStateIntergenic;

	int lastSeqLength = -1;
	int lastpos = -1;
	
	
	public GapFeaturesInterval13() {	
	}

	public int getNumFeatures() {
		return nFeatures;
	}	
	
	public String getFeatureName(int featureIndex) {
		String[] names = new String[] {"Coding frameshift", "Intron frameshift", "Intergenic frameshift", "Coding mod3 gap", "Intron mod3 gap", "Intergenic mod3 gap"};
		int raw = featureIndex - startIx;
		Assert.a(raw<nFeatures);
		String ret = names[raw];
		return ret;
	}

	transient InputSequence<? extends MultipleAlignmentColumn> lastSeq = null;
	
	public void evaluateNode(InputSequence<? extends MultipleAlignmentColumn> seq, int pos, int state, FeatureList result) {
		// Try out one of the following two lines.
		if ( (seq != lastSeq) ) {
		//if( (seq.length() != lastSeqLength)  || (pos < lastpos) ) {
			log.debug("Performing precomputations for seq of length " + seq.length() + " at position " + pos);
			performPrecomputations(seq.getX(0).getMultipleAlignment());
			lastSeqLength = seq.length();
			lastpos = pos;
			lastSeq = seq;
		}

		if (isStateCoding[state]     && frameshifter[pos])          { result.addFeature(startIx+0, 1.0); }
		if (isStateIntronic[state]   && frameshifter[pos])          { result.addFeature(startIx+1, 1.0); }
		if (isStateIntergenic[state] && frameshifter[pos])          { result.addFeature(startIx+2, 1.0); }
		if (isStateCoding[state]     && gapboundary[pos])           { result.addFeature(startIx+3, 1.0); }
		if (isStateIntronic[state]   && gapboundary[pos])           { result.addFeature(startIx+4, 1.0); }
		if (isStateIntergenic[state] && gapboundary[pos])           { result.addFeature(startIx+5, 1.0); }
	}

	private void performPrecomputations(MultipleAlignmentInputSequence seq) {
		// In this method is where you need to set the boolean vectors frameshifter[i] and gapboundary[i] using the multiple alignment seq.

		// Maybe a little inefficient with memory allocation but hopefully not too much
		if (seq.length() > frameshifter.length) {
			frameshifter = new Boolean[seq.length()];
			gapboundary  = new Boolean[seq.length()];
		}
		
		for (int j=0; j<seq.length(); j++) {
			frameshifter[j] = false;
			gapboundary[j] = false;
		}
		
		int numSpecies = seq.numSpecies();
		int consensusLength = seq.getConsensusLength();
		
		for (int spec = 0; spec<numSpecies; spec++) {
			boolean inGap = false;
			int conGapStart = -1;
			for (int cpos = 1; cpos< consensusLength; cpos++ ) {
				if ( (!inGap) && (h.hashable(seq.characterInPaddedAlignment(cpos-1,spec))) && (!h.hashable(seq.characterInPaddedAlignment(cpos,spec))) ) {
					inGap = true;
					conGapStart = cpos;
				}
				if ( (inGap) && (!h.hashable(seq.characterInPaddedAlignment(cpos-1,spec))) && (h.hashable(seq.characterInPaddedAlignment(cpos,spec))) ) {
					inGap = false;
					int conGapEnd = cpos-1;
					int gaplen = conGapEnd - conGapStart + 1;
					if (gaplen <=60) {
						if ( (gaplen%3) == 0) {
							gapboundary[seq.con2refLeft(conGapStart)] = true;
							gapboundary[seq.con2refRight(conGapEnd)] = true;
						} else {
							frameshifter[seq.con2refLeft(conGapStart)] = true;
							frameshifter[seq.con2refRight(conGapEnd)] = true;							
						}
						
					}
					
				}
					
			}
			
		}
		
	}

	public void train(int startingIndex, ModelManager modelInfo, List<? extends TrainingSequence<? extends MultipleAlignmentColumn>> data) {

		startIx = startingIndex;
		model = modelInfo;
		nStates = model.getNumStates();

		maxSeqLength = 0;			
		for(TrainingSequence<? extends MultipleAlignmentColumn> seq : data) {
			int len = seq.length();
			if (len > maxSeqLength) { maxSeqLength = len; }
		}
		frameshifter = new Boolean[maxSeqLength];
		gapboundary  = new Boolean[maxSeqLength];
				
		isStateCoding = new Boolean[nStates];       for (int j=0; j<nStates; j++) { isStateCoding[j] = false; }
		isStateCoding[1] = true;
		isStateCoding[2] = true;
		isStateCoding[3] = true;
		isStateCoding[7] = true;
		isStateCoding[8] = true;
		isStateCoding[9] = true;		

		isStateIntronic = new Boolean[nStates];     for (int j=0; j<nStates; j++) { isStateIntronic[j] = false; }
		isStateIntronic[4] = true;
		isStateIntronic[5] = true;
		isStateIntronic[6] = true;
		isStateIntronic[10] = true;
		isStateIntronic[11] = true;
		isStateIntronic[12] = true;

		isStateIntergenic = new Boolean[nStates];   for (int j=0; j<nStates; j++) { isStateIntergenic[j] = false; }
		isStateIntergenic[0] = true;
		
	}
	@Override
	public CacheStrategySpec getCacheStrategy() {
		return new CacheStrategySpec(CacheStrategy.SPARSE);
	}
	
}

