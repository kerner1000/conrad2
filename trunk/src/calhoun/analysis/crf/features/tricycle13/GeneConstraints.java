package calhoun.analysis.crf.features.tricycle13;

import java.util.List;

import calhoun.analysis.crf.AbstractFeatureManager;
import calhoun.analysis.crf.FeatureList;
import calhoun.analysis.crf.FeatureManagerEdge;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.analysis.crf.io.TrainingSequence;

/** lmplements basic constraints on gene calls.
 * <ol>
 * <li>Transition from intergenic to start must occur at ATG
 * <li>Splice sites must be canonical GT/AG or GC/AG
 * <li>Transition from exon to stop must be followed by a stop codon
 * </ol>
 */

// NOTE: I don't think these constraints look at frame, so I recommend using a version of this class adapted to the specific model yu're considering.  JPV 20060629

public class GeneConstraints extends AbstractFeatureManager<Character> implements FeatureManagerEdge<Character> {

	private static final long serialVersionUID = -3045999729941973608L;

	enum Constraint {NONE, START, DONOR, ACCEPTOR, STOP, START_MINUS, DONOR_MINUS, ACCEPTOR_MINUS, STOP_MINUS, CODING, CODING_MINUS};
	Constraint[] constraints;
	int numStates;
	
	public String getFeatureName(int featureIndex) {
		return "Gene constraints";
	}

	/** This is a constraint class, so we don't return features */
	public int getNumFeatures() {
		return 0;
	}

	/** Set up the matrix
	 * Depends on states starting with the words 'intergenic, intron, and exon'.  Also depends on the negative strand states ending in m.
	 */
	public void train(int startingIndex, ModelManager modelInfo, List<? extends TrainingSequence<? extends Character>> data) {
		numStates = modelInfo.getNumStates();
		constraints = new Constraint[numStates * numStates];
		for(int i=0; i<numStates; ++i) {
			for(int j=0; j<numStates; ++j) {
				String lastState = modelInfo.getStateName(i);
				String currentState = modelInfo.getStateName(j);
				Constraint c = Constraint.NONE;
				if(lastState.startsWith("intergenic") && currentState.startsWith("exon")) {
					c = currentState.endsWith("m") ? Constraint.STOP_MINUS : Constraint.START;
				}
				else if(lastState.startsWith("exon") && currentState.startsWith("intron")) {
					c = currentState.endsWith("m") ? Constraint.ACCEPTOR_MINUS : Constraint.DONOR;
				}
				else if(lastState.startsWith("intron") && currentState.startsWith("exon")) {
					c = currentState.endsWith("m") ? Constraint.DONOR_MINUS : Constraint.ACCEPTOR;
				}
				else if(lastState.startsWith("exon") && currentState.startsWith("intergenic")) {
					c = lastState.endsWith("m") ? Constraint.START_MINUS : Constraint.STOP;
				}
				else if(lastState.equals("exon3") && currentState.equals("exon1")) {
					c = Constraint.CODING;
				}
				else if(lastState.equals("exon1m") && currentState.equals("exon3m")) {
					c = Constraint.CODING_MINUS;
				}
				constraints[i*numStates + j] = c;
			}
		}
	}
	
	public void evaluateEdge(InputSequence<? extends Character> seq, int pos, int prevState, int state, FeatureList result) {
		boolean valid = true;

		switch(constraints[prevState*numStates + state]) {
			case START:
				valid = startConstraintPlus(seq, pos);
				break;
			case DONOR:
				valid = donorConstraintPlus(seq, pos);
				break;
			case ACCEPTOR:
				valid = acceptorConstraintPlus(seq, pos);
				break;
			case STOP:
				valid = stopConstraintPlus(seq, pos);
				break;
			case START_MINUS:
				valid = startConstraintMinus(seq, pos);
				break;
			case DONOR_MINUS:
				valid = donorConstraintMinus(seq, pos);
				break;
			case ACCEPTOR_MINUS:
				valid = acceptorConstraintMinus(seq, pos);
				break;
			case STOP_MINUS:
				valid = stopConstraintMinus(seq, pos);
				break;
			case CODING:
				valid = !stopConstraintPlus(seq, pos);
				break;
			case CODING_MINUS:
				valid = !stopConstraintMinus(seq, pos);
				break;
		}
		
		if(valid == false)
			result.invalidate();
	}
	
	boolean startConstraintPlus(InputSequence<? extends Character> seq, int pos) {
		return (seq.length() > pos + 2) && seq.getX(pos) == 'A' && seq.getX(pos+1) == 'T' && seq.getX(pos+2) == 'G';
	}

	boolean startConstraintMinus(InputSequence<? extends Character> seq, int pos) {
		return (pos >= 3) && seq.getX(pos-3) == 'C' && seq.getX(pos-2) == 'A' && seq.getX(pos-1) == 'T';
	}

	boolean donorConstraintPlus(InputSequence<? extends Character> seq, int pos) {
		return (seq.length() > pos + 1) && seq.getX(pos) == 'G' && (seq.getX(pos+1) == 'T' || seq.getX(pos+1) == 'C');
	}

	/**  CCCGTCCCAGCCC 
	 *   GGGCAGGGTCGGG
	 */
	boolean donorConstraintMinus(InputSequence<? extends Character> seq, int pos) {
		return (pos >= 2) && (seq.getX(pos-2) == 'A' || seq.getX(pos-2) == 'G') && seq.getX(pos-1) == 'C';
	}

	boolean acceptorConstraintPlus(InputSequence<? extends Character> seq, int pos) {
		return (pos > 1) && seq.getX(pos-2) == 'A' && seq.getX(pos-1) == 'G';
	}

	boolean acceptorConstraintMinus(InputSequence<? extends Character> seq, int pos) {
		return (seq.length() > pos + 1) && seq.getX(pos) == 'C' && seq.getX(pos+1) == 'T';
	}

	boolean stopConstraintPlus(InputSequence<? extends Character> seq, int pos) {
		if(seq.length() > pos + 2 && seq.getX(pos) == 'T') {
			return (seq.getX(pos+1) == 'A' && (seq.getX(pos+2) == 'G' || seq.getX(pos+2) == 'A'))
					|| (seq.getX(pos+1) == 'G' && seq.getX(pos+2) == 'A');
		}
		return false;
	}

	boolean stopConstraintMinus(InputSequence<? extends Character> seq, int pos) {
		if(pos >= 3 && seq.getX(pos-1) == 'A') {
			return (seq.getX(pos-2) == 'T' && (seq.getX(pos-3) == 'C' || seq.getX(pos-3) == 'T'))
					|| (seq.getX(pos-2) == 'C' && seq.getX(pos-3) == 'T');
		}
		return false;
	}
}
