package calhoun.analysis.crf.io;

import java.util.Collection;
import java.util.HashMap;

import calhoun.util.Assert;

/** represents an input sequence that also contains a sequence of hidden states.
 * @param <A> the elements which make up the input sequence which is part of this training sequence 
 */
public class TrainingSequence<A> implements InputSequence<A> {
	private static final long serialVersionUID = -443081006327395045L;
	InputSequence<? extends A> x;
	int[] y;
	HashMap<String,TrainingSequence> componentMap;
	
	/** constructs a training sequence using just the input sequence.  The hidden sequence needs to be added later
	 * using setY.
	 * @param xArg the input sequence the hidden states correspond to
	 */
	public TrainingSequence(InputSequence<? extends A> xArg) {
		this.x = xArg;
	}
	
	/** constructs a training sequence using these input and hidden sequences
	 * @param xArg the input sequence the hidden states correspond to
	 * @param yArg the state indices of the values of the hidden states
	 */
	public TrainingSequence(InputSequence<? extends A> xArg, int[] yArg) {
		this.x = xArg;
		this.y = yArg;
		Assert.a(y.length == x.length(), "Lengths differ between input and training sequences.  Hidden = " + y.length  + "   Observed = " + x.length());
	}
	
	/** gets the underlying input sequence
	 * @return the underlying input sequence
	 */
	public InputSequence<? extends A> getInputSequence() {
		return x;
	}

	/** gets the vector of hidden state indices
	 * @return an array of ints containing the index of the hidden state at each position in the sequence */
	public int[] getY() {
		return y;
	}
	
	/** sets the vector of hidden state indices.  Must be the same length as the previsouly specified inputs
	 * @param hiddenStates an array of ints containing the index of the hidden state at each position in the sequence */
	public void setY(int[] hiddenStates) {
		Assert.a(hiddenStates.length == x.length(), "Lengths differ between input and training sequences.  Hidden = " + hiddenStates.length  + "   Observed = " + x.length());
		y = hiddenStates;
	}
	
	/** gets the hidden state index at a particular position
	 * @param x1 the 0-based position at which to get the hidden state index
	 * @return the hidden state index at this position
	 */
	public int getY(int x1) {
		return y[x1];
	}
	
	/** sets the hidden state index at a particular position
	 * @param x the 0-based position at which to get the hidden state index
	 * @param z the vlaue of the hidden state index to set for this position
	 */
	public void setY(int x,int z) {
		y[x] = z;
	}

	/** gets the value of the underlying input sequence at a particular position
	 * @param ix a zero-based index into the input sequence
	 * @return the value of the input sequence at the position
	 */
	public A getX(int ix) {
		return x.getX(ix);
	}

	/** gets the length of the input and training sequences
	 * @return the length of the training sequence
	 */
	public int length() {
		return x == null ? 0 : x.length();
	}
	
	public InputSequence<?> getComponent(String name) {
		return getTrainingComponent(name);
	}

	public Collection<String> listComponents() {
		return x.listComponents();
	}

	/** returns a new TrainingSequence created by taking a single component of the input sequence and
	 * pairing it with the hidden states for this Training Sequence.  The input sequence must be a composite.
	 * @param name the name of the component to extract from the composite input sequence
	 * @return a new TrainingSequence created from teh extracted input sequence component
	 */
	public TrainingSequence getTrainingComponent(String name) {
		if(componentMap == null) {
			componentMap = new HashMap<String, TrainingSequence>();
		}
		TrainingSequence ret = componentMap.get(name);
		if(ret == null) {
			ret = new TrainingSequence(x.getComponent(name), y);
			componentMap.put(name, ret);
		}
		return ret;
	}

	public TrainingSequence<A> subSequence(int start, int end) {		
		Assert.a(start >= 1);
		Assert.a(end <= this.length());
		Assert.a(start <= end);
		int[] newdata = new int[end-start+1];
		for (int j=0; j<(end-start+1); j++) {
			newdata[j] = y[j+start-1];
		}
	
		TrainingSequence<A> TS = new TrainingSequence<A>(this.x.subSequence(start,end),newdata);	
		return TS;
	}

	@Override
	public String toString() {
		return x+ " + training.";
	}
}
