package calhoun.analysis.crf.features.supporting;

import java.io.Serializable;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.io.InputSequence;
import calhoun.seq.KmerHasher;
import calhoun.util.Assert;


public class PWMLookup implements Serializable {
	private static final long serialVersionUID = 5353716408013134581L;

	private static final Log log = LogFactory.getLog(PWMLookup.class);

	final KmerHasher.CharacterHash hashForward = KmerHasher.ACGTN;
	final KmerHasher.CharacterHash hashReverse = KmerHasher.ACGTNcomp;		

	boolean finalized = false;
	
	final int mult = 4;
	
	int left;
	int right;
	int span;
	
	double[] lookupTable;

	public PWMLookup(int lookLeft, int lookRight, double pseudoCount) {
		left  = lookLeft;
		right = lookRight;
		span  = left + right;
		
		Assert.a( left  >= 0);  // Perhaps all you really need is (lookLeft+lookRight) >= 0, but leave that for another day.
		Assert.a( right >= 0);
		Assert.a( span  < 30);  // No particular need for this, constraint can be relaxed, just want to be alerted if this happens for now

		int tableSize = mult*span;
	
		lookupTable = new double[tableSize];
		for (int i=0; i<tableSize; i++) {
			lookupTable[i] = pseudoCount;
		}
	}

	
	public void increment(InputSequence<? extends Character> seq, int pos, boolean isPlus) {
		Assert.a(!finalized);
		
		if (isPlus) {
			if (pos < left) { return; }
			if (pos + right > seq.length()) { return; }
			for (int j = pos - left; j<pos+right; j++) {
				int h = hashForward.hash(seq.getX(j));
				if (h<4) {
					lookupTable[mult*(j-pos+left) + h] += 1.0;
				}
			}
		} else {
			if (pos < right) { return; }
			if (pos + left > seq.length()) { return; }			
			for (int j = pos + left - 1; j>=pos - right; j--) {
				int h = hashReverse.hash(seq.getX(j));
				if (h<4) {
					lookupTable[mult*(pos+left-1-j) + h] += 1.0;
				}
			}
		}
	}

	public void completeCounts() {
		Assert.a(!finalized);
		log.debug("finalizing a PWMlookup, span="+span + "    mult=" + mult);
		
		for (int i=0; i<span; i++) {
			// System.out.println("finalizing i=" + i);
			double sum = 0;
			for (int j=mult*i; j<mult*(i+1); j++) {
				sum += lookupTable[j];
			}
			for (int j=mult*i; j<mult*(i+1); j++) {
				lookupTable[j] = Math.log(lookupTable[j]/sum);
				Assert.a(lookupTable[j]<=0);
			}
		}
		
		finalized = true;
	}
	
	public double lookup(InputSequence<? extends Character> seq, int pos, boolean isPlus) {
		Assert.a(finalized);
		
		double ret = 0.0;
		
		if (isPlus) {
			if (pos < left) { return 0.0; }
			if (pos + right > seq.length()) { return 0.0; }
			for (int j = pos - left; j<pos+right; j++) {
				int h = hashForward.hash(seq.getX(j));
				if (h<4) {
					ret += lookupTable[mult*(j-pos+left) + h];
				}
			}
		} else {
			if (pos < right) { return 0.0; }
			if (pos + left > seq.length()) { return 0.0; }			
			for (int j = pos + left - 1; j>=pos - right; j--) {
				int h = hashReverse.hash(seq.getX(j));
				if (h<4) {
					ret += lookupTable[mult*(pos+left-1-j) + h];
				}
			}
		}
		Assert.a(ret<=0);
		return ret;
	}

}