package calhoun.seq;

import java.io.Serializable;

import calhoun.util.Assert;

public class SimpleFastaSequence implements FastaSequence, Serializable {
	private static final long serialVersionUID = 3541355741077645447L;
	String header;
	private String sequence;
	private byte [] quality = null;

	public SimpleFastaSequence() {
		
	}
	
	public SimpleFastaSequence(String header, String sequence) {
		this.header = header;
		this.sequence = sequence;
	}
	
	public void setHeader(String s) {
		header = s;
	}
	
	public void setSequence(String s) {
		if(quality != null && s != null) {
			Integer sl = new Integer(s.length());
			Integer ql = new Integer(quality.length);
			Assert.a(sl.equals(ql), "Quality (",ql,") and sequence (",sl,") for ",header==null?"":header," are not the same size.");
		}
		sequence = s;
	}
	
	public String getHeader() {
		return header;
	}

	public boolean hasSequence() {
		return sequence != null;
	}

	public boolean hasQuality() {
		return quality != null;
	}

	public String getSequence() {
		return sequence;
	}

	public String getSequence(int start, int stop) {
		Assert.a(sequence != null, "No sequence available for ", header);
		return sequence.substring(start-1, stop);
	}

	public int getLength() {
		if(sequence != null) {
			return sequence.length();
		}
		else {
			Assert.a(quality != null, "No sequence or quality available for ", header);
			return quality.length;
		}
	}
	
	public byte[] getQuality() {
		return quality;
	}
	
	public void setQuality(byte [] q) {
		if(q != null && sequence != null) {
			Assert.a(sequence.length() == q.length, header, ": Quality (",new Integer(q.length),") and sequence (",new Integer(sequence.length()),") for ",header," are not the same size.");
		}
		quality = q;
	}
}
