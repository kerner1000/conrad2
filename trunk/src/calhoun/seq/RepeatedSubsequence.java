package calhoun.seq;

import calhoun.util.Assert;
import calhoun.util.PrimeUtil;

public class RepeatedSubsequence {

	/** Test if infinite repetitions of a and b can be aligned */
	public static boolean isRollingMatch(String a, String b) {
		int len = a.length();
		if(len != b.length()) {
			return false;
		}
		if(a.equals(b)) {
			return true;
		}
		StringBuffer sb = new StringBuffer(a);
		for(int i = 1; i<len; ++i) {
			for(int j=0; j<len;++j) {
				sb.setCharAt(j, a.charAt((i+j)%len));
			}
			if(sb.toString().equals(b)) {
				return true;
			}
		}
		return false;
	}
	
	/** Returns the largest repeated subsequence in a given string. */
	public static String calc(String seq) {
		
		// Pos is the end of the shortest repeat found so far.  starts as being the whole length of the sequence.
		int pos = seq.length();
	
		// Prime our indexes with the original sequence
		int maxPrime = seq.length()/2+1;
		
		// Get list of all possible prime factors
		int[] primes = PrimeUtil.primesToAtLeastN(maxPrime);
		
		// Iterate up the list of primes (start at 2)
		for(int i = 1; primes[i] <= maxPrime; ++i) {
			// While the current repeat length is divisible by the current prime
			while(pos%primes[i] == 0) {
				boolean repeat = checkRepeat(seq, pos/primes[i], primes[i]);
				if(repeat) {
					// If it is identical, choose the first section and attempt to split that
					pos = pos/primes[i];
				}
				else {
					// If it is not, try the next prime
					break;
				}
			}
		}

		String ret = seq.substring(0, pos);
		// Also check for all the same if we got this far
		if(checkRepeat(ret, 1, ret.length())) {
			ret = ret.substring(0,1);
		}
		
		// Verify results.
		Assert.a(checkRepeat(seq, ret.length(), seq.length()/ret.length()));
		return ret;
	}

	/** Checks to see if the given sequence begins with n repeats, each of length l. */
	public static boolean checkRepeat(String s, int l, int n) {
		for(int i = 0; i<l; ++i) {
			char c = s.charAt(i);
			for(int j = 1; j<n; ++j) {
				if (c != s.charAt(i+j*l))
					return false;
			}
		}
		return true;
	}
}
