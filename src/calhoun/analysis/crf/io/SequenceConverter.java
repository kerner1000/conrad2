package calhoun.analysis.crf.io;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import calhoun.util.Assert;

/** used for converting between different encodings of gene structure in a hidden sequence.  Mostly for legacy use
*/
public class SequenceConverter {

	private static HashMap<String, Integer> map = new HashMap<String, Integer>();

	
	public static ArrayList<ArrayList<Integer>> stateVector2StateLengths(List<? extends TrainingSequence<?>> data, int nStates) {

		ArrayList<ArrayList<Integer>> durations = new ArrayList<ArrayList<Integer>>();
		
		for (int j=0; j<nStates; j++) {
			durations.add(new ArrayList<Integer>());
		}
		
		for (TrainingSequence<?> seq : data) {
			if (seq.length()==0) { continue; }
			int oldState = seq.getY(0);
			int intervalStart = 0;
			for (int pos=1; pos<seq.length(); pos++) {
				int newState = seq.getY(pos);
				if (newState != oldState) {
					durations.get(oldState).add(pos-intervalStart);
					intervalStart = pos;
					oldState = newState;
				}
			}
			durations.get(oldState).add(seq.length()-intervalStart);
		}
		
		return durations;
	}
	
	
	public static int[] convertSeqFromInterval13ToInterval29(int[] states) {
		int swap, preswap;
		int ctr = 0;
		for(int i=0; i<states.length-1; i++) {
			swap = 0; preswap = 0;
			// ig-e
			if (states[i]==0 && states[i+1]>=1 && states[i+1]<=3) {
				preswap = 13;
			// e-ig
			} else if (states[i]>=1 && states[i]<=3 && states[i+1]==0) {
				swap = 14;
			// e-i_I
			} else if (states[i]>=1 && states[i]<=3 && states[i+1]>=4 && states[i+1]<=6) {
				swap = 15 + (states[i+1] - 4);
			// i-e_E
			} else if (states[i]>=4 && states[i]<=6 && states[i+1]>=1 && states[i+1]<=3) {
				preswap = 18 + (states[i+1] - 1);
			// ig-em
			} else if (states[i]==0 && states[i+1]>=7 && states[i+1]<=9) {
				preswap = 21;
			// em-ig
			} else if (states[i]>=7 && states[i]<=9 && states[i+1]==0) {
				swap = 22;								
			// em-i_Im
			} else if (states[i]>=7 && states[i]<=9 && states[i+1]>=10 && states[i+1]<=12) {
				swap = 23 + (states[i+1] - 10);
			// im-e_Em
			} else if (states[i]>=10 && states[i]<=12 && states[i+1]>=7 && states[i+1]<=9) {
				preswap = 26 + (states[i+1] - 7);
			}
			if (swap != 0) {
				if (i != states.length-2 && states[i+2] == states[i+1]) {
					states[i+1] = swap;
					states[i+2] = swap;
				} else if (i == states.length-2) {
					states[i+1] = swap;
				}

				ctr++;
			} else if (preswap != 0) {
				if (i >= 1 && states[i] == states[i-1]) {
					states[i] = preswap;
					states[i-1] = preswap;
				} else if (i == 0) {
					states[i] = preswap;
				}
			}
		}
		//log.debug("made " + ctr + " changes, where length is " + states.length);
		return states;
	}
	
	public static int[] convertSeqFromInterval29ToInterval13(int[] seq) {
		int len = seq.length;
		Assert.a(len>=1);
		int prevInterval29y = -1;
		for (int pos = 0; pos<len; pos++) {
			int interval29y = seq[pos];
			int interval13y = -1;
			if (interval29y == 13 || interval29y == 14 || interval29y == 21 || interval29y == 22) {
				interval13y = 0;
			} else if (interval29y >= 15 && interval29y <= 17) {
				interval13y = (interval29y - 15) + 4;
			} else if (interval29y >= 23 && interval29y <= 25) {
				interval13y = (interval29y - 23) + 10;
			} else if (interval29y >= 18 && interval29y <= 20) {
				Assert.a(prevInterval29y >= 4 && prevInterval29y <= 6);
				interval13y = prevInterval29y;
			} else if (interval29y >= 26 && interval29y <= 28) {
				if (prevInterval29y == -1) {
					// XXX: hack
					prevInterval29y = interval29y - 26 + 10;
				}
				Assert.a(prevInterval29y >= 10 && prevInterval29y <= 12, "prev is " + prevInterval29y);
				interval13y = prevInterval29y;
			} else {
				interval13y = interval29y;
			}
			Assert.a(interval13y != -1);
			Assert.a(interval13y < 13);
			seq[pos] = interval13y;
			if (interval29y < 18 || (interval29y > 20 && interval29y < 26)) {
				prevInterval29y = interval29y;
			}
		}	
		return seq;
	}
	
	public static int[] convertSeqFromInterval29ToInterval13Wrong(int[] states) {
		int swap;
		for(int i=0; i<states.length; i++) {
			swap = states[i];
			if (states[i] == 13 || states[i] == 14 || states[i] == 21 || states[i] == 22) {
				swap = 0;
			} else if (states[i] >= 15 && states[i] <= 17) {
				swap = states[i] - 15 + 4;
			// Likely wrong
			} else if (states[i] >= 18 && states[i] <= 20) {
				swap = states[i] - 18 + 4;
			} else if (states[i] >= 23 && states[i] <= 25) {
				swap = states[i] - 23 + 10;
			// Also likely wrong
			} else if (states[i] >= 26 && states[i] <= 28) {
				swap = states[i] - 26 + 10;
			}
			
			Assert.a(swap >= 0 && swap <= 12);
			states[i] = swap;
		}
		return states;
	}
	
	
	public static void convertSeqFromTricycle13ToInterval13(TrainingSequence<Character> seq) {
		int len = seq.length();
		Assert.a(len>=1);
		
		for (int pos = 0; pos<len; pos++) {
			int tricycle13y = seq.getY(pos);
			int interval13y = posTricycle2interval13(pos,tricycle13y);
			seq.setY(pos,interval13y);
		}
	}

	public static String convertSeqFromTricycle13ToInterval13(String seq2) {
		char[] seq = seq2.toCharArray();
		int len = seq.length;
		Assert.a(len>=1);
		
		for (int pos = 0; pos<len; pos++) {
			int tricycle13y = char2integer13(seq[pos]);
			int interval13y = posTricycle2interval13(pos,tricycle13y);
			seq[pos] = integer132char(interval13y);
		}
		return new String(seq);
	}
	
	private static char integer132char(int i) {
		if (i<10) {
			return (char) ('0'+i);
		} else if (i<36) {
			return (char) ('A'+(i-10));
		} else if (i<62) {
			return (char) ('a'+(i-36));
		}
		Assert.a(false);
		return 0;
	}


	private static int char2integer13(char x) {		
		int temp = x - '0';
		if ( (temp<0) || (temp>9)) {
			temp = x - 'A' + 10;			
			if ( (temp<10) || (temp>35)) {
				temp = x - 'a' + 36;
				Assert.a( (temp>=36) && (temp<62), "Offending character was '" + x);
			}
		}
		Assert.a(temp<13,"temp = " + temp + "  and x = " + x);
		return temp;
	}
	
	private static int posTricycle2interval13(int pos, int tricycle13y) {
		int interval13y = 0;
		switch(tricycle13y) {
		case 0:
			// intergenic
			interval13y = 0;
			break;
		case 1: 
		case 2:
		case 3:
			// exon plus strand
			interval13y = ((( pos - (tricycle13y - 1) ) %3 +3) %3) + 1;
			break;
		case 4:
			// intron plus strand
			interval13y = 6;
			break;
		case 5:
			interval13y = 5;
			break;
		case 6:
			interval13y = 4;
			break;
		case 7:
		case 8:
		case 9:
			// exon minus strand
			interval13y = (pos + (tricycle13y - 7) + 1) % 3 + 7;
			break;
		case 10:
			// intron minus strand
			interval13y = 12;
			break;
		case 11:
			interval13y = 11;
			break;
		case 12:
			interval13y = 10;
			break;
		default:
			Assert.a(false);
		}
		return interval13y;
	}


	// NOTE: it just so happens that the conversion Interval13ToTricycle13 is it's own inverse,
	// but this is just a coincidence and best programming prcatice to write it out twice.
	public static void convertSeqFromInterval13ToTricycle13(TrainingSequence<Character> seq) {
		int len = seq.length();
		Assert.a(len>=1);
		
		for (int pos = 0; pos<len; pos++) {
			int interval13y = seq.getY(pos);
			int tricycle13y = posInterval2tricycle13(pos,interval13y);
			seq.setY(pos,tricycle13y);
		}	
	}	
	
	public static String convertSeqFromInterval13ToTricycle13(String seq2) {
		char[] seq = seq2.toCharArray();
		int len = seq.length;
		Assert.a(len>=1);
		
		for (int pos = 0; pos<len; pos++) {
			int interval13y = char2integer13(seq[pos]);
			int tricycle13y = posInterval2tricycle13(pos,interval13y);
			seq[pos] = integer132char(tricycle13y);
		}	
		return new String(seq);
	}		
	
	
	private static int posInterval2tricycle13(int pos, int interval13y) {
		int tricycle13y = 0;
		switch(interval13y) {
		case 0:
			// intergenic
			tricycle13y = 0;
			break;
		case 1: 
		case 2:
		case 3:
			// exon plus strand
			tricycle13y = ((( pos - (interval13y - 1) ) %3 +3) %3) + 1;
			break;
		case 4:
			// intron plus strand
			tricycle13y = 6;
			break;
		case 5:
			tricycle13y = 5;
			break;
		case 6:
			tricycle13y = 4;
			break;
		case 7:
		case 8:
		case 9:
			// exon minus strand
			tricycle13y = (((-pos + (interval13y - 7) + 2) %3 +3) %3) + 7;
			break;
		case 10:
			// intron minus strand
			tricycle13y = 12;
			break;
		case 11:
			tricycle13y = 11;
			break;
		case 12:
			tricycle13y = 10;
			break;
		default:
			Assert.a(false);
		}
		return tricycle13y;
	}






	// Given a hidden sequence that has 13 states and values [0, 12], converts that
	// sequence to a 39 state model with values [0, 38].
	public static void convertSeqFrom13To39(TrainingSequence<Character> seq)
	{
		setStateMap();
		if (seq.length() < 2)	return;
		
		ArrayList<SeqPair> states = new ArrayList<SeqPair>();
		int i, state39, total, k, seqIdx, curIdx;
		int seqLen = seq.length();
		int startElement = seq.getY(0);
		int startIndex   = 0;

		int prevElement = seq.getY(0);
		int curElement  = seq.getY(1);
		int nextElement;
		int prevState   = -1;
		
		for (i=2; i<seqLen; i++)
		{		
			curIdx = i-1;
			
			nextElement = seq.getY(i);
			Assert.a(curElement>=0 && curElement <=12, "invalid character in hidden sequence, '", curElement, "'");				
			
			if (!sameState(startElement, curElement))
			{
				// End the previous state
				state39 = getState39(startElement, prevElement, curElement, prevState);
				states.add(new SeqPair(state39, (curIdx-1)-startIndex + 1));
				prevState = state39;
				
				// Start the current state
				startElement = curElement;
				startIndex   = curIdx;
			}
			prevElement = curElement;
			curElement  = nextElement;
		}

		// Add last state
		state39 = getState39(startElement, prevElement, -1, prevState);
		states.add(new SeqPair(state39, (i-1)-startIndex+1));

		// Verify sequence lengths will be the same
		total = 0;
		for (i=0; i<states.size(); i++)
			total += states.get(i).length;
		Assert.a(total == seqLen, "Sum of state lengths = " + total + ", Sequence Length = " + seqLen);
		
		// Set the values in the sequence to be in 39 state model.
		seqIdx = 0;
		for (i=0; i<states.size(); i++)
		{
			for (k=0; k<states.get(i).length; k++)
			{
				seq.setY(seqIdx, states.get(i).state);
				seqIdx++;
			}
		}
	}

	
	// Given a hidden sequence that has 39 states and values [0, 38], converts that
	// sequence to a 13 state model with values [0, 12].
	public static void convertSeqFrom39To13(TrainingSequence<Character> seq)
	{
		int len = seq.length();
		if (len < 1)	{return;}
		
		int[] y = new int[len];
		
		for (int j=0; j<len; j++) {
			y[j] = seq.getY(j); }
		
		convertSeqFrom39To13(y);
		
		for (int j=0; j<len; j++) {
			seq.setY(j,y[j]); }
		
	}
	
	// Given a hidden sequence that has 39 states and values [0, 38], converts that
	// sequence to a 13 state model with values [0, 12].
	public static void convertSeqFrom39To13(int[] seq)
	{
		if (seq.length < 1)	{return;}
		
		int seqLen = seq.length;
		int cur, i, exonPhase;
		boolean inExon = false;
		exonPhase = -1;
		
		cur = seq[0];
		for (i=1; i<seqLen; i++)
		{		
			cur = seq[i];
			Assert.a(cur>=0 && cur <=38, "invalid character in hidden sequence, '", cur, "'");		

			if (cur==0)					// INTERGENIC, do nothing, 0 in both models
			{
				inExon = false;
			}
			else if (isIntron39(cur))	// INTRON
			{
				seq[i] = convertIntron39To13(cur);
				inExon = false;
			}
			else 								// EXON
			{
				if (inExon)		// already been here, just keep cycling through 1, 2, 3, or 9, 8, 7
				{
					seq[i] = exonPhase;
					exonPhase = incrementExonPhase(exonPhase);
				}
				else			// first time we're entering an exon
				{
					inExon = true;
					exonPhase = convertExon39To13(cur);
					seq[i] = exonPhase;
					exonPhase = incrementExonPhase(exonPhase);
				}
			}
		}
	}	
	
	
	//
	// SUPPORTING FUNCTIONS FOR 13 -> 39 CONVERSION
	//
	
	// Given a start and end state in the 13 state model, returns the state in the 39 state model.
	private static int getState39(int start, int end, int next, int prevState)
	{
		int state39 = -1;
		
		if      (start ==  0 && end ==  0)	state39 = map.get("NTG").intValue();
		else if (start ==  6 && end ==  6)	state39 = map.get("I0p").intValue();
		else if (start ==  4 && end ==  4)	state39 = map.get("I1p").intValue();
		else if (start ==  5 && end ==  5)	state39 = map.get("I2p").intValue();
		else if (start == 12 && end == 12)	state39 = map.get("I0m").intValue();
		else if (start == 10 && end == 10)	state39 = map.get("I1m").intValue();
		else if (start == 11 && end == 11)	state39 = map.get("I2m").intValue();

		else if (start ==  1 && end ==  3)	
		{
			if (prevState == map.get("NTG").intValue() && next == 0)	state39 = map.get("ENNp").intValue();
			else if (prevState == map.get("NTG").intValue())			state39 = map.get("EN0p").intValue();
			else if (next == 0)											state39 = map.get("E0Np").intValue();
			else														state39 = map.get("E00p").intValue();
		}
		else if (start ==  1 && end ==  1)	
		{
			if (prevState == map.get("NTG").intValue())		state39 = map.get("EN1p").intValue();
			else											state39 = map.get("E01p").intValue();
		}
		else if (start ==  1 && end ==  2)	
		{
			if (prevState == map.get("NTG").intValue())		state39 = map.get("EN2p").intValue();
			else											state39 = map.get("E02p").intValue();
		}
		else if (start ==  2 && end ==  3)	
		{
			if (next == 0)									state39 = map.get("E1Np").intValue();
			else 											state39 = map.get("E10p").intValue();
		}
		else if (start ==  3 && end ==  3)	
		{
			if (next == 0)									state39 = map.get("E2Np").intValue();
			else 											state39 = map.get("E20p").intValue();
		}
		else if (start ==  2 && end ==  1)	state39 = map.get("E11p").intValue();
		else if (start ==  2 && end ==  2)	state39 = map.get("E12p").intValue();
		else if (start ==  3 && end ==  1)	state39 = map.get("E21p").intValue();
		else if (start ==  3 && end ==  2)	state39 = map.get("E22p").intValue();
		
		else if (start ==  9 && end ==  7)	
		{
			if (prevState == map.get("NTG").intValue() && next == 0)	state39 = map.get("ENNm").intValue();
			else if (prevState == map.get("NTG").intValue())			state39 = map.get("E0Nm").intValue();
			else if (next == 0)											state39 = map.get("EN0m").intValue();
			else 														state39 = map.get("E00m").intValue();
		}
		else if (start ==  7 && end ==  7)	
		{
			if (next == 0)									state39 = map.get("EN1m").intValue();
			else 											state39 = map.get("E01m").intValue();
		}
		else if (start ==  8 && end ==  7)	
		{
			if (next == 0)									state39 = map.get("EN2m").intValue();
			else 											state39 = map.get("E02m").intValue();
		}
		else if (start ==  9 && end ==  8)	
		{
			if (prevState == map.get("NTG").intValue())		state39 = map.get("E1Nm").intValue();
			else											state39 = map.get("E10m").intValue();
		}
		else if (start ==  9 && end ==  9)	
		{
			if (prevState == map.get("NTG").intValue())		state39 = map.get("E2Nm").intValue();
			else 											state39 = map.get("E20m").intValue();
		}
		else if (start ==  7 && end ==  8)	state39 = map.get("E11m").intValue();
		else if (start ==  8 && end ==  8)	state39 = map.get("E12m").intValue();
		else if (start ==  7 && end ==  9)	state39 = map.get("E21m").intValue();
		else if (start ==  8 && end ==  9)	state39 = map.get("E22m").intValue();
		
		if (state39 == -1)	Assert.a(false, "start = " + start + "   end = " + end);
		return state39;
	}
	
	// Given a two states in the 13 state model, returns true if they are the same state,
	// else returns false.  I.e. 1,2,3 are all positive exons and considered the same state.
	private static boolean sameState(int state1, int state2)
	{
		if (state1 == state2)
			return true;
		else if (isPlusExon(state1) && isPlusExon(state2))
			return true;
		else if (isMinusExon(state1) && isMinusExon(state2))
			return true;
		return false;
	}
	
	// Given a state in the 13 state model, returns true if it is an exon on the plus strand.
	private static boolean isPlusExon(int state)
	{
		if (state == 1 || state == 2 || state == 3)
			return true;
		return false;
	}
	
	// Given a state in the 13 state model, returns true if it is an exon on the minus strand.
	private static boolean isMinusExon(int state)
	{
		if (state == 9 || state == 8 || state == 7)
			return true;
		return false;
	}
	
	//
	// SUPPORTING FUNCTIONS FOR 39 -> 13 CONVERSION
	//
	
	// Given an INTRON state in the 39 model state, returns the intron state in the 13 model state.
	private static int convertIntron39To13(int element)
	{	
		if (element == map.get("I1p").intValue()) 	return 4;	// intron1
		if (element == map.get("I2p").intValue()) 	return 5;	// intron2
		if (element == map.get("I0p").intValue()) 	return 6;	// intron3
		if (element == map.get("I1m").intValue()) 	return 10;	// intron1m
		if (element == map.get("I2m").intValue()) 	return 11;	// intron2m
		if (element == map.get("I0m").intValue()) 	return 12;	// intron3m
		
		Assert.a(false);
		return -1;
	}
	
	// Given an EXON state in the 39 model state, returns the intron state in the 13 model state.
	private static int convertExon39To13(int element)
	{
		int exonPhase = -1;
		
		if (getStrand39(element) == +1)		// plus strand
		{
			if      (element == map.get("ENNp").intValue())		exonPhase = 1;	// exon1
			else if (element == map.get("EN0p").intValue())		exonPhase = 1;	// exon1
			else if (element == map.get("EN1p").intValue())		exonPhase = 1;	// exon1
			else if (element == map.get("EN2p").intValue())		exonPhase = 1;	// exon1
			else if (element == map.get("E00p").intValue())		exonPhase = 1;	// exon1
			else if (element == map.get("E01p").intValue())		exonPhase = 1;	// exon1
			else if (element == map.get("E02p").intValue())		exonPhase = 1;	// exon1
			else if (element == map.get("E10p").intValue())		exonPhase = 2;	// exon2
			else if (element == map.get("E11p").intValue())		exonPhase = 2;	// exon2
			else if (element == map.get("E12p").intValue())		exonPhase = 2;	// exon2
			else if (element == map.get("E20p").intValue())		exonPhase = 3;	// exon3
			else if (element == map.get("E21p").intValue())		exonPhase = 3;	// exon3
			else if (element == map.get("E22p").intValue())		exonPhase = 3;	// exon3
			else if (element == map.get("E0Np").intValue())		exonPhase = 1;	// exon1
			else if (element == map.get("E1Np").intValue())		exonPhase = 2;	// exon2
			else if (element == map.get("E2Np").intValue())		exonPhase = 3;	// exon3
		}
		else if (getStrand39(element) == -1)	// minus strand
		{
			if      (element == map.get("ENNm").intValue())		exonPhase = 9;	// exon3m
			else if (element == map.get("EN0m").intValue())		exonPhase = 9;	// exon3m
			else if (element == map.get("EN1m").intValue())		exonPhase = 7;	// exon1m
			else if (element == map.get("EN2m").intValue())		exonPhase = 8;	// exon2m
			else if (element == map.get("E00m").intValue())		exonPhase = 9;	// exon3m
			else if (element == map.get("E01m").intValue())		exonPhase = 7;	// exon1m
			else if (element == map.get("E02m").intValue())		exonPhase = 8;	// exon2m
			else if (element == map.get("E10m").intValue())		exonPhase = 9;	// exon3m
			else if (element == map.get("E11m").intValue())		exonPhase = 7;	// exon1m
			else if (element == map.get("E12m").intValue())		exonPhase = 8;	// exon2m
			else if (element == map.get("E20m").intValue())		exonPhase = 9;	// exon3m
			else if (element == map.get("E21m").intValue())		exonPhase = 7;	// exon1m
			else if (element == map.get("E22m").intValue())		exonPhase = 8;	// exon2m
			else if (element == map.get("E0Nm").intValue())		exonPhase = 9;	// exon3m
			else if (element == map.get("E1Nm").intValue())		exonPhase = 9;	// exon3m
			else if (element == map.get("E2Nm").intValue())		exonPhase = 9;	// exon3m
		}
		if (exonPhase == -1)	Assert.a(false);
		return exonPhase;
	}
	
	// Given an exon in the 39 state model, returns the strand (+1 or -1)
	private static int getStrand39(int element)
	{
		if (7 <= element && element <= 22)
			return (+1);
		else if (23 <= element && element <= 38)
			return (-1);
		
		Assert.a(false);
		return 0;
	}
	
	// Given an element in the 39 model state, returns true if this element
	//   is an intron, else returns false
	private static boolean isIntron39(int element)
	{
		if (1 <= element && element <= 6)
			return true;
		else
			return false;
	}
	
	private static int incrementExonPhase(int phase)
	{
		if (phase == 1 || phase == 2 || phase == 3)
		{
			phase++;
			if (phase == 4)	phase = 1;
		}
		else if (phase == 9 || phase == 8 || phase == 7)
		{
			phase--;
			if (phase == 6)	phase = 9;
		}
		else
		{
			Assert.a(false);
		}
		return phase;
	}
	
	// 
	// CONVERT A HIDDEN SEQUENCE TO A GFF FILE
	// 
	
	// This function converts a 13 state model hidden sequence to a GFF file.  
	// Used for debugging only.
	public static void writeHiddenSequenceGFF(TrainingSequence<Character> refStates, String filename) throws IOException
	{
		int i, ref, exonStart, exonEnd, geneNum;
		boolean inPlusGene = false;
		boolean inMinusGene = false;
		geneNum = 0;
		exonStart = exonEnd = geneNum = 0;
		Writer fout = new BufferedWriter(new FileWriter(filename));	

		for (i=0; i<refStates.length(); i++)
		{
			ref = refStates.getY(i);
			
			if (ref == 1 || ref == 2 || ref == 3)
			{
				if (!inPlusGene)
				{
					exonStart = i+1;
					inPlusGene = true;
				}
			}
			else if (ref == 7 || ref == 8 || ref == 9)
			{
				if (!inMinusGene)
				{
					exonStart = i+1;
					inMinusGene = true;
				}
			}
			else if (inPlusGene && (ref == 4 || ref == 5 || ref == 6) )
			{
				exonEnd = i;
				fout.write("XXX\tghmm\tENNp\t" + exonStart + "\t" + exonEnd + "\t0\t+\t0\tgene_" + geneNum + "\n");
				inPlusGene = false;
			}
			else if (inMinusGene && (ref == 10 || ref == 11 || ref == 12) )
			{
				exonEnd = i;
				fout.write("XXX\tghmm\tENNm\t" + exonStart + "\t" + exonEnd + "\t0\t-\t0\tgene_" + geneNum + "\n");
				inMinusGene = false;
			}
			else
			{
				if (inPlusGene || inMinusGene)	// was in gene at last nucleotide, but not in any more
				{
					exonEnd = i;
					
					if (inPlusGene)
						fout.write("XXX\tghmm\tENNp\t" + exonStart + "\t" + exonEnd + "\t0\t+\t0\tgene_" + geneNum + "\n");
					else 
						fout.write("XXX\tghmm\tENNm\t" + exonStart + "\t" + exonEnd + "\t0\t-\t0\tgene_" + geneNum + "\n");
					inPlusGene = false;
					inMinusGene = false;
					geneNum++;
				}
			}
			
		}
		fout.close();
	}
	
	// This function converts a 39 state model hidden sequence to a GFF file.  
	// Used for debugging only.
	public void writeHiddenSequence39GFF(TrainingSequence<Character> refStates, String filename) throws IOException
	{
		int i, ref, exonStart, exonEnd, geneNum;
		boolean inPlusGene = false;
		boolean inMinusGene = false;
		geneNum = 0;
		exonStart = exonEnd = geneNum = 0;
		Writer fout = new BufferedWriter(new FileWriter(filename));	
		for (i=0; i<refStates.length(); i++)
		{
			ref = refStates.getY(i);
			
			if (isPlusExon39(ref))	// plus exon
			{
				if (!inPlusGene)
				{
					exonStart = i+1;
					inPlusGene = true;
				}
			}
			else if (isMinusExon39(ref))	// minus exon
			{
				if (!inMinusGene)
				{
					exonStart = i+1;
					inMinusGene = true;
				}
			}
			else if (inPlusGene && isPlusIntron39(ref) )	// plus intron
			{
				exonEnd = i;
				fout.write("XXX\tghmm\t" + stateIdxToString(refStates.getY(exonStart)) + "\t" + exonStart + "\t" + exonEnd + "\t0\t+\t0\tgene_" + geneNum + "\n");
				inPlusGene = false;
			}
			else if (inMinusGene && isMinusIntron39(ref) )	// minus intron
			{
				exonEnd = i;
				fout.write("XXX\tghmm\t" + stateIdxToString(refStates.getY(exonStart)) + "\t" + exonStart + "\t" + exonEnd + "\t0\t-\t0\tgene_" + geneNum + "\n");
				inMinusGene = false;
			}
			else
			{
				if (inPlusGene || inMinusGene)	// was in gene at last nucleotide, but not in any more
				{
					exonEnd = i;
					
					if (inPlusGene)
						fout.write("XXX\tghmm\t" + stateIdxToString(refStates.getY(exonStart)) + "\t" + exonStart + "\t" + exonEnd + "\t0\t+\t0\tgene_" + geneNum + "\n");
					else 
						fout.write("XXX\tghmm\t" + stateIdxToString(refStates.getY(exonStart)) + "\t" + exonStart + "\t" + exonEnd + "\t0\t-\t0\tgene_" + geneNum + "\n");
					inPlusGene = false;
					inMinusGene = false;
					geneNum++;
				}
			}
			
		}
		fout.close();	
	}
	
	// Returns true if state is a exon in the 39 state model, else returns false
	private boolean isPlusExon39(int state)
	{
		if (7 <= state && state <= 22)
			return true;
		return false;
	}
	private boolean isMinusExon39(int state)
	{
		if (23 <= state && state <= 38)
			return true;
		return false;
	}
	private boolean isPlusIntron39(int state)
	{
		if (1 <= state && state <= 3)
			return true;
		return false;
	}
	private boolean isMinusIntron39(int state)
	{
		if (4 <= state && state <= 6)
			return true;
		return false;
	}	
	
	// 
	// COMMON FUNCTIONS AND STRUCTS
	//
	
	// NOTE:  This info is copied from GHMM.  It likely should not be in two places.
	private static void setStateMap()
	{
		int statenum = 0;
		
		map.put("NTG", new Integer(statenum));	statenum++;
		map.put("I0p", new Integer(statenum));	statenum++;
		map.put("I1p", new Integer(statenum));	statenum++;
		map.put("I2p", new Integer(statenum));	statenum++;
		map.put("I0m", new Integer(statenum));	statenum++;
		map.put("I1m", new Integer(statenum));	statenum++;
		map.put("I2m", new Integer(statenum));	statenum++;
		map.put("ENNp", new Integer(statenum));	statenum++;
		map.put("EN0p", new Integer(statenum));	statenum++;
		map.put("EN1p", new Integer(statenum));	statenum++;
		map.put("EN2p", new Integer(statenum));	statenum++;
		map.put("E00p", new Integer(statenum));	statenum++;
		map.put("E01p", new Integer(statenum));	statenum++;
		map.put("E02p", new Integer(statenum));	statenum++;
		map.put("E10p", new Integer(statenum));	statenum++;
		map.put("E11p", new Integer(statenum));	statenum++;
		map.put("E12p", new Integer(statenum));	statenum++;
		map.put("E20p", new Integer(statenum));	statenum++;
		map.put("E21p", new Integer(statenum));	statenum++;
		map.put("E22p", new Integer(statenum));	statenum++;
		map.put("E0Np", new Integer(statenum));	statenum++;
		map.put("E1Np", new Integer(statenum));	statenum++;
		map.put("E2Np", new Integer(statenum));	statenum++;
		map.put("ENNm", new Integer(statenum));	statenum++;
		map.put("EN0m", new Integer(statenum));	statenum++;
		map.put("EN1m", new Integer(statenum));	statenum++;
		map.put("EN2m", new Integer(statenum));	statenum++;
		map.put("E00m", new Integer(statenum));	statenum++;
		map.put("E01m", new Integer(statenum));	statenum++;
		map.put("E02m", new Integer(statenum));	statenum++;
		map.put("E10m", new Integer(statenum));	statenum++;
		map.put("E11m", new Integer(statenum));	statenum++;
		map.put("E12m", new Integer(statenum));	statenum++;
		map.put("E20m", new Integer(statenum));	statenum++;
		map.put("E21m", new Integer(statenum));	statenum++;
		map.put("E22m", new Integer(statenum));	statenum++;
		map.put("E0Nm", new Integer(statenum));	statenum++;
		map.put("E1Nm", new Integer(statenum));	statenum++;
		map.put("E2Nm", new Integer(statenum));	statenum++;
	}	
	
	private String stateIdxToString(int state)
	{
		switch(state)
		{
		case 0: return "NTG";
		case 1: return "I0p";
		case 2: return "I1p";
		case 3: return "I2p";
		case 4: return "I0m";
		case 5: return "I1m";
		case 6: return "I2m"; 
		case 7: return "ENNp";
		case 8: return "EN0p";
		case 9: return "EN1p";
		case 10: return "EN2p";
		case 11: return "E00p";
		case 12: return "E01p";
		case 13: return "E02p";
		case 14: return "E10p";
		case 15: return "E11p";
		case 16: return "E12p";
		case 17: return "E20p";
		case 18: return "E21p";
		case 19: return "E22p";
		case 20: return "E0Np";
		case 21: return "E1Np";
		case 22: return "E2Np";
		case 23: return "ENNm";
		case 24: return "EN0m";
		case 25: return "EN1m";
		case 26: return "EN2m";
		case 27: return "E00m";
		case 28: return "E01m";
		case 29: return "E02m";
		case 30: return "E10m";
		case 31: return "E11m";
		case 32: return "E12m";
		case 33: return "E20m";
		case 34: return "E21m";
		case 35: return "E22m";
		case 36: return "E0Nm";
		case 37: return "E1Nm";
		case 38: return "E2Nm";
		}
		return "XXX";
	}
	
	private static class SeqPair
	{
		public int state;
		public int length;
		public SeqPair(int st, int len) {state=st; length=len; }
	}


}
