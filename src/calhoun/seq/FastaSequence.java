package calhoun.seq;

/** 
 * Represents a generic sequence in a fasta file.  FastaIterator and FastaReader
 * have different classes that they store the result in, but these methods are 
 * common.  Unifies the Fasta classes so you can more easily intermix them. 
 */
public interface FastaSequence {
	/** The number of sequence characters */
	public int getLength();
	
	/** Gets the fasta header (with leading and trailing spaces trimmed) */
	public String getHeader();
	
	/** Returns the entire sequence.  */
	public String getSequence();
	
	/** Returns a portions of the sequence.
	 * @param start 1-based index of the first character to return.  Can be from 1 to the length of the sequence.
	 * @param stop 1-based index of the last charater to return.  Can be from 1 to the length of the sequence.
	 */
	public String getSequence(int start, int stop);

	/** Returns the quality scores for this sequence.  May be null if no quality is present. */ 
	public byte[] getQuality();
}
