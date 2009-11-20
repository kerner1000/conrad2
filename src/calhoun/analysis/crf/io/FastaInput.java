package calhoun.analysis.crf.io;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import calhoun.seq.FastaIterator;
import calhoun.seq.FastaSequence;
import calhoun.seq.FastaWriter;
import calhoun.util.Assert;

public class FastaInput implements InputComponentIO {
	private static final long serialVersionUID = 760405914814389112L;

	String header;
	String sequence;
	Integer numSeqCap = Integer.MAX_VALUE;

	public List<String> getComponentNames() {
		List<String> ret = new ArrayList();
		ret.add(header);
		ret.add(sequence);
		return ret;
	}

	public void readInputSequences(String location, List<Map<String, InputSequence<?>>> inputs) throws IOException {
		FastaIterator it = new FastaIterator(location);
		Assert.a(inputs.size() == 0, "Fasta input always has to go first in the list.");
		int i = 0;
		if (numSeqCap == null) {
			numSeqCap = Integer.MAX_VALUE;
		}
		while(it.hasNext() && i < numSeqCap.intValue()) {
			FastaSequence seq = (FastaSequence) it.next();
			Map<String, InputSequence<?>> input = new HashMap();
			input.put(header, new NameInputSequence(seq.getHeader()));
			input.put(sequence, new InputSequenceCharacter(seq.getSequence().toUpperCase()));
			inputs.add(input);
			i++;
		}
	}

	public void writeInputSequences(String location, List<? extends Map<String, ? extends InputSequence<?>>> inputComponents) throws IOException {
		FastaWriter w = new FastaWriter(location, false);
		for(Map<String, ? extends InputSequence<?>> input : inputComponents) {
			NameInputSequence nameSeq = (NameInputSequence) input.get(header);
			InputSequenceCharacter seqSeq = (InputSequenceCharacter) input.get(sequence);
			w.writeSeq(nameSeq.getName(), seqSeq.getString());
		}
		
		w.close();
	}

	/**
	 * @return Returns the header.
	 */
	public String getHeader() {
		return header;
	}
	/**
	 * @param header The header to set.
	 */
	public void setHeader(String header) {
		this.header = header;
	}
	/**
	 * @return Returns the sequence.
	 */
	public String getSequence() {
		return sequence;
	}
	/**
	 * @param sequence The sequence to set.
	 */
	public void setSequence(String sequence) {
		this.sequence = sequence;
	}
	
	public Integer getNumSeqCap() {
		return numSeqCap;
	}
	
	public void setNumSeqCap(Integer numSeqCap) {
		this.numSeqCap = numSeqCap; 
	}
}
