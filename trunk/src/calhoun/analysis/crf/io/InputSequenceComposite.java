package calhoun.analysis.crf.io;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import calhoun.util.Assert;

/** a composite input sequence made up of individual components.  This is useful when putting together different
 * features that may take different inputs.  Each component is given a name within the overall composite and can
 * be reference by that name.  All sequences must have the same length.  The value at any position in the sequence
 * is a map that relates component names to their values at that position.  */
public class InputSequenceComposite implements InputSequence<Map<String, Object>> {

	Map<String, InputSequence<?>> components = new HashMap<String, InputSequence<?>>();
	int length = -1;

	/** default constructor.  Components should be added with addComponents. */
	public InputSequenceComposite() {
	}
	
	/** constructor that initializes the sequence with a set of components. 
	 * @param data the initial input sequences to use as the components */
	public InputSequenceComposite(Map<String, InputSequence<?>> data) {
		addComponents(data);
	}
	
	/** adds a new component input sequence to this composite.  Each new component must have a unique name and all 
	 * components must have the same length. */
	public void addComponent(String name, InputSequence<?> seq) {
		Assert.a(!components.containsKey(name), "Component already exists in composite input: ", name);
		if(length == -1) {
			length = seq.length();
		}
		Assert.a(seq.length() == -1 || length == seq.length(), "Composite sequence has length ", length, " but attempted to add component named ", name, " of length ", seq.length());
		components.put(name, seq);
	}
	
	/** adds all input sequences in the given map as components. */
	public void addComponents(Map<String, InputSequence<?>> data) {
		for(Map.Entry<String, InputSequence<?>> entry : data.entrySet()) {
			addComponent(entry.getKey(), entry.getValue());
		}
	}

	/** returns all of components of this input sequence.
	 * @return a map that maps component names to input sequences
	 */
	public Map<String, InputSequence<?>> getComponents() {
		return components;
	}

	public Map<String, Object> getX(int ix) {
		Map<String, Object> ret = new HashMap<String, Object>();
		for(Map.Entry<String, InputSequence<?>> entry : components.entrySet()) {
			ret.put(entry.getKey(), entry.getValue().getX(ix));
		}
		return ret;
	}

	public int length() {
		Assert.a(components.size() > 0, "Length of composite sequence was called before any components were added.");
		return length;
	}

	public InputSequence<?> getComponent(String name) {
		return components.get(name);
	}

	public Collection<String> listComponents() {
		return components.keySet();
	}

	public InputSequenceComposite subSequence(int start, int end) {
		Assert.a(start >= 1);
		Assert.a(end <= length());
		Assert.a(start <= end);		
		
		InputSequenceComposite ret = new InputSequenceComposite();
		for(Map.Entry<String, InputSequence<?>> entry : components.entrySet()) {
			ret.addComponent(entry.getKey(), entry.getValue().subSequence(start,end));
		}
		
		return ret;
	}
}
