package calhoun.analysis.crf.features.supporting;

import java.io.Serializable;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.io.InputSequence;
import calhoun.seq.KmerHasher;
import calhoun.util.Assert;


public class LogProbLookup implements Serializable {
	private static final long serialVersionUID = -9195647924401633963L;
	private static final Log log = LogFactory.getLog(LogProbLookup.class);
	final KmerHasher.CharacterHash hashForward = KmerHasher.ACGTother;
	final KmerHasher.CharacterHash hashReverse = KmerHasher.ACGTotherRC;		
	boolean finalized;
	
	final int mult = 4;
	int maxLookBack;
	
	double[] lookupTable;
	int lookupTableSize;
	
	public LogProbLookup(int lookBack, double pseudoCount) {
		Assert.a(lookBack >= 0);
		Assert.a(lookBack < 10);
		this.maxLookBack = lookBack;
		finalized = false;
		
		lookupTableSize = 1;
		for (int i=0; i<=lookBack; i++) {
			lookupTableSize *= mult;
		}
		lookupTable = new double[lookupTableSize];
		
		for (int i=0; i<lookupTableSize; i++) {
			lookupTable[i] = pseudoCount;
		}
	}

	private boolean isHistory(InputSequence<? extends Character> seq, int pos) {

		for (int j=pos-maxLookBack; j<=pos+maxLookBack; j++) {
			if (hashForward.hash(seq.getX(j))==4) { return false; }
			// Above is identical to checking hashReverse
			// If there are N's within history window in either dircetion, want to ignore this position
		}
		return true;
	}
	
	private int getInd(InputSequence<? extends Character> seq, int pos, boolean isPlus) {
		int ind = 0;
		if (isPlus) {
			if (pos < maxLookBack) { return -1; }
			for (int j=pos-maxLookBack; j<=pos; j++) {
				int h = hashForward.hash( (char) seq.getX(j));
				//int h = hashForward.hash('A');
				if (h<4) {
					ind *= mult;
					ind += h;
				} else {
					return -1;
				}
			}
		} else {
			if (pos + maxLookBack >= seq.length()) { return -1; }
			for (int j=pos+maxLookBack; j>=pos; j--) {
				int h = hashReverse.hash( (char) seq.getX(j));
				if (h<4) {
					ind *= mult;
					ind += h;
				} else {
					return -1;
				}
			}
		}
		return ind;
	}
	
	public void increment(InputSequence<? extends Character> seq, int pos, boolean isPlus) {
		Assert.a(!finalized);
		int ind = getInd(seq,pos,isPlus);
		if (ind >=0) {
			lookupTable[ind] += 1.0;
		}
	}

	@Override
	public void finalize() {
		Assert.a(!finalized);
		log.debug("finalizing a LogProbLookup, lookupTablesize="+lookupTableSize + "    mult=" + mult);
		for (int i=0; i<lookupTableSize/mult; i++) {
			// System.out.println("finalizing i=" + i);
			double sum = 0;
			for (int j=mult*i; j<mult*(i+1); j++) {
				sum += lookupTable[j];
			}
			for (int j=mult*i; j<mult*(i+1); j++) {
				lookupTable[j] = Math.log(lookupTable[j]/sum);
			}
		}
		finalized = true;
	}
	
	public double lookup(InputSequence<? extends Character> seq, int pos, boolean isPlus) {
		Assert.a(finalized);	
		int ind = getInd(seq,pos,isPlus);
		if (ind >= 0) {
			return lookupTable[ind];
		}
		return 0.0;
	}

}