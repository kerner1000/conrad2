package calhoun.analysis.crf.solver.check;

import calhoun.analysis.crf.ModelManager;
import calhoun.util.Assert;
import calhoun.util.DenseBooleanMatrix2D;
import calhoun.util.DenseIntMatrix2D;

public class TransitionInfo {
	public ModelManager fm;
	public short nStates;
	public int nTransitions;
	public int nPotentials;

	/// Valid transitions
	public short[] transitionFrom;
	public short[] transitionTo;
	// A matrix containing the sparse index of the given transition.  -1 for transitions not allowed in the model
	public DenseIntMatrix2D transitionIndex;
	// An ordering of potentials that puts each node potential directly before all the edge potentials leading into that node
	public short[] orderedPotentials;
	
	public boolean ignoreSemiMarkovSelf;
	/** Fill the transitionFrom and transitionTo arrays for a fast Mi lookup */
	public TransitionInfo(ModelManager fm, boolean allPaths) {
		initTrans(fm, allPaths);
	}
	
	/** Fill the transitionFrom and transitionTo arrays for a fast Mi lookup */
	protected TransitionInfo() {
	}
	
	/** returns true if self transitions should be allowed to the given state, even if not present in the model. */
	protected boolean allowSelf(int state) {
		return false;
	}
	
	protected void initTrans(ModelManager argFm, boolean allPaths) {
		this.fm = argFm;
		nStates = (short) argFm.getNumStates();
		transitionIndex = new DenseIntMatrix2D(nStates, nStates);
		transitionIndex.assign(-1);
		DenseBooleanMatrix2D transitions = argFm.getLegalTransitions();
		if(transitions == null || allPaths) {
			transitions = new DenseBooleanMatrix2D(nStates, nStates);
			transitions.assign(true);
		}

		short count = 0;
		for(short i = 0; i<nStates; ++i) {
			// This should probably be in a derived class, but is here because it is an important check I don't have a better place for now.
			Assert.a(ignoreSemiMarkovSelf || !allowSelf(i) || !transitions.getQuick(i, i), "Self transitions are not currently allowed in the semi-Markov model.  Set the 'ignoreSemiMarkovSelfTransitions' property of the CacheProcessor to true if you want to ignore this problem.  This is useful when using the same edge and transition set for a semi-markov and regular markov model.", i);
			for(short j = 0; j<nStates; ++j) {
				if(transitions.getQuick(i, j) || (i == j && allowSelf(i)))
					count++;
			}
		}
		nTransitions = count;
		nPotentials = nStates + nTransitions;
		orderedPotentials = new short[nPotentials];
		transitionFrom = new short[nTransitions];
		transitionTo = new short[nTransitions];
		count = 0;
		int orderedCount = 0;
		for(short i = 0; i<nStates; ++i) {
			orderedPotentials[orderedCount] = i;
			orderedCount++;
			for(short j = 0; j<nStates; ++j) {
				if(transitions.getQuick(j, i) || (i == j && allowSelf(i))) {
					orderedPotentials[orderedCount] = (short) (nStates+count);
					orderedCount++;
					transitionIndex.setQuick(j, i, count);
					transitionFrom[count] = j;
					transitionTo[count] = i;
					++count;
				}
			}
		}
		Assert.a(count == nTransitions);
	}
}
