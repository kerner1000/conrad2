package calhoun.analysis.crf.io;

import java.util.Collection;

import calhoun.util.Assert;

/** input sequence that captures the 'name' of a sequence.  This is useful for display or for tracking outputs,
 * but usually won't be used by the engine.  The full sequence name is returned as the data element for every position.
 * In this respect it is not really a sequence.<p>
 * This feature has some special properties that allow sequences to be tracked as they are subsetted.
 * When subsetting a <code>NameInputSequence</code>, a string is append to the name showing the start and stop (one-based
 * inclusive) of the subsetting.  For example, if the sequence with name "CND1" is subsetting to based 201 to 1000, then
 * the new name will be "CND1:201-1000".  To avoid repetitive subsetting, the NameInput does not allow subsetting on names
 * containing a '-'.
 */
public class NameInputSequence implements InputSequence<String>  {
	String str;

	/** default constructor */
	public NameInputSequence() { }
	
	/** constructors an input sequence using this name */
	public NameInputSequence(String a) {
		str = a;
	}
	
	public String getName() {
		return str;
	}

	public String getX(int ix) {
		return str;
	}

	public int length() {
		return -1;
	}

	public InputSequence<?> getComponent(String name) {
		throw new UnsupportedOperationException();
	}

	public Collection<String> listComponents() {
		throw new UnsupportedOperationException();
	}

	public InputSequence<String> subSequence(int start, int end)  {
		Assert.a(!str.contains("-"));
		NameInputSequence n = new NameInputSequence(str + ":" + start + "-" + end);
		return n;
	}
}