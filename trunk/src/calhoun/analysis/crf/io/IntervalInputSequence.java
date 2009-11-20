package calhoun.analysis.crf.io;

import java.util.Collection;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.util.Assert;
import calhoun.util.RangeMap;

/** an input sequence made up of intervals of constant valued features.  Each position is represented
 * by an {@link IntervalPosition} object.
 */
public class IntervalInputSequence implements InputSequence<IntervalInputSequence.IntervalPosition>{
	@SuppressWarnings("unused")
	private static final Log log = LogFactory.getLog(IntervalInputSequence.class);

	RangeMap rmplus,rmminus;
	String inputName;
	int inputLength;
	
	/** constructs an input sequence from the necessary data
	 * @param rmplus range map of values on the positive strand
	 * @param rmminus range map of values on the negative strand
	 * @param inputName name of the input sequence
	 * @param inputLength length of the input sequence.  Requires since the intervals may not cover the whole sequence.
	 */
	public IntervalInputSequence(RangeMap rmplus, RangeMap rmminus, String inputName, int inputLength) {
		this.rmplus = rmplus;
		this.rmminus = rmminus;
		this.inputName = inputName;
		this.inputLength = inputLength;
	}
	
	public IntervalPosition getX(int ix) {
		IntervalPosition ret = new IntervalPosition(ix);
		return ret;
	}

	public int length() {
		return inputLength;
	}
	
	public IntervalInputSequence subSequence(int start, int end) { // 1-based inclusive
		Assert.a(end >= start);
		Assert.a(end <= inputLength);
		Assert.a(start >= 1);
		
		RangeMap newRmplus  = shiftedSubRangeMap(rmplus ,start,end);
		RangeMap newRmminus = shiftedSubRangeMap(rmminus,start,end);		

		return new IntervalInputSequence(newRmplus,newRmminus,inputName,end-start+1);
	}

	public InputSequence getComponent(String name) {
		throw new UnsupportedOperationException();
	}

	public Collection<String> listComponents() {
		throw new UnsupportedOperationException();
	}

	private RangeMap shiftedSubRangeMap(RangeMap rm, int start, int end) {
		Object[] rmplusList = rm.find(start-1,end-1).toArray();
		RangeMap newrm = new RangeMap();
		for (int i=0; i<rmplusList.length; i++) {
			IntervalRangeMapValue irmv = (IntervalRangeMapValue) rmplusList[i];
			IntervalRangeMapValue newirmv = new IntervalRangeMapValue(irmv.start-start+1, irmv.end-start+1, irmv.value);
			rm.add(newirmv.start,newirmv.end, newirmv);
		}
		
		return newrm;
	}

	/** provides interval based queries for a particular position in the input sequence for an IntervalInputSequence. 
	 */
	public class IntervalPosition {
		int pos;
		
		/** creates an IntervalPosition object for this sequence at this position. 
		 * @param pos the position on the input sequence */
		public IntervalPosition(int pos) {
			this.pos = pos;
		}

		/** returns true if an input interval exists on the positive strand at this position.  */
		public boolean queryPlus() {
			return (rmplus.find(pos,pos+1).size()>0);
		}
		
		/** returns true if an input interval exists on the negative strand at this position.  */
		public boolean queryMinus() {
			return (rmminus.find(pos,pos+1).size()>0);
		}
	}

	static public class IntervalRangeMapValue {
		public int start;
		public int end;
		public double value;
		
		public IntervalRangeMapValue(int start, int end, double value) {
			this.start = start;
			this.end = end;
			this.value = value;
		}
		
		public void insertIntoRangeMap(RangeMap RM) {
			RM.add(start,end,this);
		}
		
		public String toStringStrand(String strand) {
			String ret = "(" + start + "," + end + "," + strand + "," + value + ")";
			return ret;
		}
		
	}
}
