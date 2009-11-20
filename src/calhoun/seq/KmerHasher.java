package calhoun.seq;

import java.io.Serializable;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.util.Assert;
import calhoun.util.ErrorException;

/** Class for computing kmer hashes.  You instantiate the class with the length of the kmer to hash and the alphabet to use.  
 * The hashing functions then compute hashes for individual kmers.  This class does not stored any hashes, which are just ints. */
public class KmerHasher implements Serializable {
	private static final long serialVersionUID = -3402947063680917230L;

	@SuppressWarnings("unused")
	private static final Log log = LogFactory.getLog(KmerHasher.class);

	CharacterHash charHash;
	int length;
	int mscMult;	// Multiplier for most significant character
	int lscMult;	// Multiplier for least significant character
	
	public interface CharacterHash extends Serializable {
		public short getSize();
		public int hash(char a) throws ErrorException;
		//public boolean hashable(char a);
		public char reverse(int a);
	}
	
	/** Creates a hash from a given character hash and length.
	 */
	public KmerHasher(CharacterHash charHash, int length) {
		this.charHash = charHash;
		this.length = length;
		mscMult = (int) Math.pow(charHash.getSize(), length-1);
		lscMult = charHash.getSize();
		// log.info(String.format("Alphabet size: %d Max hash: %d Max hash(len-1): %d, length: %d", lscMult, mscMult*lscMult-1, mscMult, length));
	}

	
	public static String reverseComplement(String forward) {
		String reverse = "";
		
		CharacterHash hf = ACGTother;
		CharacterHash hr = ACGTotherRC;
		
		for (int j=forward.length()-1; j>=0; j--) {
			reverse = reverse + hr.reverse(hf.hash(forward.charAt(j)));
		}
		
		return reverse;
	}
	
	
	public int range() {
		//return (int) mscMult+lscMult-1;//Math.pow(charHash.getSize(), length);
		return (int) mscMult*lscMult;
	}
	
	/** Computes hash given a 0-based position on the string. */
	public int hash(String str, int pos) {
		Assert.a(pos>=0);
		Assert.a(pos+length <= str.length());
		int hash = 0;
		for(int i = 0; i<length; ++i) {
			hash = hash*lscMult + charHash.hash(str.charAt(pos+i));
		}
		return hash;
	}

	/** Computes a hash given a character array.  The array must be the exact size of the hash */
	public int hash(char[] chr) {
		int hash = 0;
		for(int i = 0; i<length; ++i) {
			hash = hash*lscMult + charHash.hash(chr[i]);
		}
		return hash;
	}

	/** Computes a hash given a character array.  Starts at the given index into the array */
	public int hash(char[] chr, int start) {
		int hash = 0;
		for(int i = start; i<length+start; ++i) {
			hash = hash*lscMult + charHash.hash(chr[i]);
		}
		return hash;
	}

	/** Updates an existing hash.  Drops the first character and adds in the new one to the end.
	 * hash("BCDE", 0) == shiftHash("E", hash("ABCD", 0))*/
	public int shiftHash(char chr, int hash) {
		return ((hash%mscMult)*lscMult) + charHash.hash(chr);
		
	}

	/** Updates an existing hash.  Drops the first character and adds in the new one to the end.
	 * hash("ABCD", 0) == reverseShiftHash("A", hash("BCDE", 0))*/
	public int reverseShiftHash(char chr, int hash) {
		return charHash.hash(chr)*mscMult + hash/lscMult;
	}

	/** Character hash function to use with DNA bases.  Upper and lower case get hashed to the same value.  Handles only the 4 nucleotides (ACTG).  No other characters allowed. */
	public static CharacterHash DNA = new CharacterHash() {
		private static final long serialVersionUID = -5641887174464060367L;
		final char[] BASES = new char[] {'A','C','G','T'};
		public short getSize() { return 4;}
		public int hash(char a) {
			switch(a) {
				case 'A':
				case 'a':
					return 0;
				case 'C':
				case 'c':
					return 1;
				case 'G':
				case 'g':
					return 2;
				case 'T':
				case 't':
					return 3;
				default:
					throw new ErrorException("Bad character for hashing '"+a+"'.  Only A,C,T,G,a,c,t,g are allowed.");
			}
			
		}
		public char reverse(int a) {
			return BASES[a];
		}
	};
	
	/** Character hash function to use with DNA bases which included the ambiguity code "N".  Upper and lower case get hashed to the same value.  Handles only the 4 nucleotides (ACTG).  No other characters allowed. */
	public static CharacterHash ACGTN = new CharacterHash() {
		private static final long serialVersionUID = -4366495850512839766L;
		final char[] BASES = new char[] {'A','C','G','T','N'};
		public short getSize() { return 5;}
		public int hash(char a) {
			switch(a) {
				case 'A':
				case 'a':
					return 0;
				case 'C':
				case 'c':
					return 1;
				case 'G':
				case 'g':
					return 2;
				case 'T':
				case 't':
					return 3;
				case 'N':
				case 'n':
					return 4;
				default:
					throw new ErrorException("Bad character for hashing '"+a+"'.  Only A,C,T,G,N,a,c,t,g,n are allowed.");
			}
		}
		public char reverse(int a) {
			return BASES[a];
		}
	};
	
	/** Character hash function to use with DNA bases which included the ambiguity code "N".  Upper and lower case get hashed to the same value.  Handles only the 4 nucleotides (ACTG).  No other characters allowed. */
	public static CharacterHash ACGTNcomp = new CharacterHash() {
		private static final long serialVersionUID = -5959780907260989652L;
		final char[] BASES = new char[] {'T','G','C','A','N'};
		public short getSize() { return 5;}
		public int hash(char a) {
			switch(a) {
				case 'T':
				case 't':
					return 0;
				case 'G':
				case 'g':
					return 1;
				case 'C':
				case 'c':
					return 2;
				case 'A':
				case 'a':
					return 3;
				case 'N':
				case 'n':
					return 4;
				default:
					throw new ErrorException("Bad character for hashing '"+a+"'.  Only A,C,T,G,N,a,c,t,g,n are allowed.");
			}
		}
		public char reverse(int a) {
			return BASES[a];
		}
	};
	
	
	/** Character hash function to use with DNA bases which included the ambiguity code "N".  Upper and lower case get hashed to the same value.  Handles only the 4 nucleotides (ACTG).  No other characters allowed. */
	public static CharacterHash ACGTother = new CharacterHash() {
		private static final long serialVersionUID = -3279138437137318988L;
		final char[] BASES = new char[] {'A','C','G','T','N'};
		public short getSize() { return 5;}
		public int hash(char a) {
			switch(a) {
				case 'A':
				case 'a':
					return 0;
				case 'C':
				case 'c':
					return 1;
				case 'G':
				case 'g':
					return 2;
				case 'T':
				case 't':
					return 3;
				default:
					return 4;
			}
		}
		public char reverse(int a) {
			return BASES[a];
		}
	};
	
	
	/** Character hash function to use with DNA bases which included the ambiguity code "N".  Upper and lower case get hashed to the same value.  Handles only the 4 nucleotides (ACTG).  No other characters allowed. */
	public static CharacterHash ACGTotherRC = new CharacterHash() {
		private static final long serialVersionUID = 609468521914363959L;
		/* Like DNA, excapt never throws an exception; returns 4 if not ACGTN. */
		final char[] BASES = new char[] {'T','G','C','A','N'};
		public short getSize() { return 5;}
		public int hash(char a) {
			switch(a) {
				case 'A':
				case 'a':
					return 3;
				case 'C':
				case 'c':
					return 2;
				case 'G':
				case 'g':
					return 1;
				case 'T':
				case 't':
					return 0;
				default:
					return 4;
			}
		}
		public char reverse(int a) {
			return BASES[a];
		}
	};
	
	
	/** Character hash function to use with any letters.  Upper and lower case get hashed to the same value. */
	public static CharacterHash LETTERS = new CharacterHash() { 
		private static final long serialVersionUID = -4036999685428397071L;
		public short getSize() { return 26;}
		public int hash(char a) {
			if(a >= 'a' && a <= 'z') {
				return a - 'a';
			}
			else if(a >= 'A' && a <= 'Z') {
				return a - 'A'; 
			}
			else {
				throw new ErrorException("Bad character for hashing '"+a+"'.  Only a-z and A-Z are allowed.");
			}
		}
		public char reverse(int a) {
			return (char)('a'+(char)a);
		}
	};

	public int hash(char c) {
		return charHash.hash(c);
	}

	public boolean hashable(char a)  {
		try {
			charHash.hash(a);
		} catch (ErrorException E) {
			return false;
		}
		return true;
	}
	
}
