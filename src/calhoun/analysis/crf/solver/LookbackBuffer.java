package calhoun.analysis.crf.solver;

import java.util.Arrays;

/**
 * This object holds information about previous positions during the computation of betas and expectations. This
 * allows us to quickly access data about previous positions. These objects are kept in a recycling buffer that
 * keeps one buffer for each possible lookback.
 * <p>
 * One tricky aspect of this is that the details change slightly between the forward and backwards pass.  On the forward
 * pass, the lookback contains the information in the normal way.  In the backwards pass, stable states and transitions are 
 * shifted back one base compared to the betas.
 */
public final class LookbackBuffer {
	public int pos;
	
	/// In the beta pass, the mi matrix for transitioning from pos+l to pos
	public double[] mi;

	/// For each transition to a semi-markov state, stores the probability of all segments including that transition.
	public double[] transitionProb;

	/// The weighted sum of feature values for staying in this position from the end of the sequence to this position
	public double[] stableState;
	
	/// Beta values at this position
	public double[] beta;

	/// Norm of beta values
	public int betaNorm;
	
	public LookbackBuffer(int states, int transitions) {
		mi = new double[states+transitions];
		transitionProb = new double[transitions];
		stableState = new double[states];
		beta = new double[states];
	}
	
	/** mi and stableStates are cleared as new values are entered. This fixes the others */
	public void clear()
	{
		pos = -1;
		Arrays.fill(beta, 0.0);
		betaNorm = Integer.MIN_VALUE;
		Arrays.fill(transitionProb, 0.0);
	}
}
