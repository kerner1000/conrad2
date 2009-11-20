package calhoun.analysis.crf.io;

import java.io.BufferedReader;
import java.io.IOException;

public class Interval21HiddenSequenceTranslator extends IntInput {
	private static final long serialVersionUID = 4413724139445660883L;
	
	@Override
	public int[] readSequence(BufferedReader r) throws IOException {
		int[] states = super.readSequence(r);
		// Weird boundary conditions here
		int swap;
		int ctr = 0;
		// HEY!
		for(int i=0; i<states.length-1; i++) {
			swap = 0;
			// ig-e
			if (states[i]==0 && states[i+1]>=1 && states[i+1]<=3) {
				swap = 13;
			// e-i
			} else if (states[i]>=1 && states[i]<=3 && states[i+1]>=4 && states[i+1]<=6) {
				swap = 14;
			// i-e
			} else if (states[i]>=4 && states[i]<=6 && states[i+1]>=1 && states[i+1]<=3) {
				swap = 15;
			// e-ig
			} else if (states[i]>=1 && states[i]<=3 && states[i+1]==0) {
				swap = 16;
			// ig-em
			} else if (states[i]==0 && states[i+1]>=7 && states[i+1]<=9) {
				swap = 17;
			// em-im
			} else if (states[i]>=7 && states[i]<=9 && states[i+1]>=10 && states[i+1]<=12) {
				swap = 18;
			// im-em
			} else if (states[i]>=10 && states[i]<=12 && states[i+1]>=7 && states[i+1]<=9) {
				swap = 19;
			// em-ig
			} else if (states[i]>=7 && states[i]<=9 && states[i+1]==0) {
				swap = 20;				
			}
			if (swap != 0) {
				if (i != states.length-1 && states[i+2] == states[i+1]) {
					states[i+1] = swap;
					states[i+2] = swap;
				} else if (i == states.length-1) {
					states[i+1] = swap;
				}

				ctr++;
			}
		}
		//System.out.println("made " + ctr + " changes, where length is " + states.length);
		return states;
	}
	
}
