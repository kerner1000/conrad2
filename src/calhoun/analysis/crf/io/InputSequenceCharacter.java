package calhoun.analysis.crf.io;

import java.util.Collection;

import calhoun.util.Assert;

/** an input sequence where the elements of the sequence are string characters. */
public class InputSequenceCharacter implements InputSequence<Character> {
	String str;
	
	/** default constructor */
	public InputSequenceCharacter() { }
	
	/** creates an input sequence using data from this string
	 * @param a string containing characters of the input sequence */
	public InputSequenceCharacter(String a) {
		str = a;
	}
	
	/** returns the input sequence as a string 
	 * @return a string of the characters */
	public String getString() {
		return str;
	}
	
	public Character getX(int ix) {
		return str.charAt(ix);
	}

	public int length() {
		return str == null ? 0 : str.length();
	}

	public InputSequence<?> getComponent(String name) {
		throw new UnsupportedOperationException();
	}

	public Collection<String> listComponents() {
		throw new UnsupportedOperationException();
	}

	public InputSequence<Character> subSequence(int start, int end) {
		Assert.a(start >= 1);
		Assert.a(end <= this.length());
		Assert.a(start <= end);
		
		InputSequenceCharacter S = new InputSequenceCharacter(str.substring(start-1,end));
		
		return S;
	}
	
	@Override
	public String toString() {
		return str.substring(0, Math.min(str.length(), 15));
	}
}
