package calhoun.analysis.crf.features.interval29;

import java.util.List;

import calhoun.analysis.crf.AbstractFeatureManager;
import calhoun.analysis.crf.CacheStrategySpec;
import calhoun.analysis.crf.FeatureList;
import calhoun.analysis.crf.FeatureManagerEdge;
import calhoun.analysis.crf.FeatureManagerNode;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.CacheStrategySpec.CacheStrategy;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.util.Assert;

/** Implements basic constraints on gene calls.
 * 1) Intergenic - start must occur at ATG
 * 2) Splice sites must be canonical GT/AG or GC/AG
 * 3) Exon-stop must be followed by a start codon
 */
public class GeneConstraintsInterval29  extends AbstractFeatureManager<Character> implements FeatureManagerEdge<Character>, FeatureManagerNode<Character> {
	private static final long serialVersionUID = 3041359216265032511L;
	
	public String getFeatureName(int featureIndex) {
		return "Gene constraints for the model Interval29";
	}

	/** This is a constraint class, so we don't return features */
	public int getNumFeatures() {
		return 0;
	}

	public void train(int startingIndex, ModelManager modelInfo, List<? extends TrainingSequence<? extends Character>> data) {
		Interval29Tools.verify(modelInfo);
	}
	
	public void evaluateEdge(InputSequence<? extends Character> seq, int pos, int prevState, int state, FeatureList result) {
		boolean valid = true;
		
		int eind,iind;
		int eind1, eind2, iind1, iind2;
		switch(Interval29Tools.edgeConstraints[prevState*Interval29Tools.numStates + state]) {
		case NONE:
			break;
		case NEVER:
			Assert.a(false);
			break;
		case PSTART:
			// ig-e to exonE
			eind = Interval29Tools.check012(state-1);
			if ((pos-eind)%3 != 0) { valid = false; break; } 
			valid = startConstraintPlus(seq, pos);
			break;
		case PDON:
			// exonE to e-iI 
			// k = e-i
			iind = Interval29Tools.check012(state-15);
			eind = Interval29Tools.check012(prevState-1);
			if ((pos-eind+iind)%3 != 0) { valid = false; break; } 
			valid = pos < 0 || donorConstraintPlus(seq, pos);
			break;
		case PACC:
			// intronI to i-eE
			// k = e-i
			eind = Interval29Tools.check012(state-18);
			iind = Interval29Tools.check012(prevState-4);
			if ((pos+2-eind+iind)%3 != 0) { valid = false; break; } 
			valid = (pos+2 >= seq.length()) || acceptorConstraintPlus(seq, pos+2);
			break;
		case PSTOP:
			// exonE to e-ig
			eind = Interval29Tools.check012(prevState-1);
			if ((pos-eind)%3 != 0) { valid = false; break; } 
			valid = stopEdgeConstraintPlus(seq, pos);
			break; //done to here
		case MSTART:
			// exonEm to em-ig
			eind = Interval29Tools.check012(prevState-7);
			if ((pos-eind)%3 != 0) { valid = false; break; } 
			valid = startConstraintMinus(seq, pos);
			break;
		case MDON:
			// intronIm to im-eEm
			// k = e+i
			eind = Interval29Tools.check012(state-26);
			iind = Interval29Tools.check012(prevState-10);
			if ((pos+2-eind-iind)%3 != 0) { valid = false; break; } 
			valid = (pos+2 >= seq.length()) || donorConstraintMinus(seq, pos+2);
			break;
		case MACC:
			// exonEm to e-iIm
			// k = e+i
			iind = Interval29Tools.check012(state-23);
			eind = Interval29Tools.check012(prevState-7);		
			if ((pos-eind-iind)%3 != 0) { valid = false; break; } 			
			valid = pos < 0 || acceptorConstraintMinus(seq, pos);
			break;
		case MSTOP:
			// ig-em to exonEm
			eind = Interval29Tools.check012(state - 7);
			if ((pos-eind)%3 != 0) { valid = false; break; } 
			valid = stopEdgeConstraintMinus(seq, pos);
			break;
		case PCODE: // redundant with node invalidation below
			eind = Interval29Tools.check012(state-1);
			if ( (pos-eind)%3 == 2) {
				valid = !stopNodeConstraintPlus(seq, pos);
			}
			break;
		case MCODE: // redundant with node evaluation below
			eind = Interval29Tools.check012(state-7);
			if ( (pos-eind)%3==0) {
				valid = !stopNodeConstraintMinus(seq, pos);
			}
			break;
		case PKEEPE:
			// i-e_j to e_j
			eind1 = Interval29Tools.check012(prevState - 18);
			eind2 = Interval29Tools.check012(state - 1);
			if (eind1 != eind2) {
				valid = false;
			} else {
				valid = acceptorConstraintPlus(seq, pos);
			}
			//valid = false;
			break;
		case PKEEPI:
			// e-i_j to i_j
			iind1 = Interval29Tools.check012(prevState - 15);
			iind2 = Interval29Tools.check012(state - 4);
			if (iind1 != iind2) {
				valid = false;
			} else {
				valid = pos-2 < 0 || donorConstraintPlus(seq, pos-2);
			}		
			//valid = false;
			break;
		case MKEEPE:
			// im-e_jm to e_jm
			eind1 = Interval29Tools.check012(prevState - 26);
			eind2 = Interval29Tools.check012(state - 7);
			if (eind1 != eind2) {
				valid = false;
			} else {
				valid = donorConstraintMinus(seq, pos);
			}
			//valid = false;
			break;
		case MKEEPI:
			// em-i_jm to i_jm
			iind1 = Interval29Tools.check012(prevState - 23);
			iind2 = Interval29Tools.check012(state - 10);
			if (iind1 != iind2) {
				valid = false;
			} else {
				valid = pos-2 < 0 || acceptorConstraintMinus(seq, pos-2);
			}	
			//valid = false;
			break;
		case PSTOPPED:
			valid = pos-2 < 0 || stopEdgeConstraintPlus(seq, pos-2);
			break;
		case MSTARTED:
			valid = pos-2 < 0 || startConstraintMinus(seq, pos-2);
			break;
		case PWILLSTART:
			valid = pos+2 >= seq.length() || startConstraintPlus(seq, pos+2);
			break;
		case MWILLSTOP:
			valid = pos+2 >= seq.length() || stopEdgeConstraintMinus(seq, pos+2);
			break;
		default:
			Assert.a(false);
		}
		
		if(valid == false)
			result.invalidate();
	}

	public void evaluateNode(InputSequence<? extends Character> seq, int pos, int state, FeatureList result) {
		boolean valid = true;
	
		int eind;
		
		switch(Interval29Tools.nodeConstraints[state]) {
		case NONE:
			break;
		case NEVER:
			Assert.a(false);
			break;
		case PCODE:
			eind = Interval29Tools.check012(state-1);
			if ( (pos-eind)%3 == 2) {
				valid = !stopNodeConstraintPlus(seq, pos);
			}
			break;
		case MCODE:
			eind = Interval29Tools.check012(state-7);
			if ( (pos-eind)%3==0) {
				valid = !stopNodeConstraintMinus(seq, pos);
			}
			break;
		default:
			Assert.a(false);
		}	
		if(valid == false)
			result.invalidate();
	}
	
	
	private boolean startConstraintPlus(InputSequence<? extends Character> seq, int pos) {
		return (seq.length() > pos + 2) && seq.getX(pos) == 'A' && seq.getX(pos+1) == 'T' && seq.getX(pos+2) == 'G';
	}

	private boolean startConstraintMinus(InputSequence<? extends Character> seq, int pos) {
		return (pos >= 3) && seq.getX(pos-3) == 'C' && seq.getX(pos-2) == 'A' && seq.getX(pos-1) == 'T';
	}

	private boolean donorConstraintPlus(InputSequence<? extends Character> seq, int pos) {
		return (seq.length() > pos + 1) && seq.getX(pos) == 'G' && (seq.getX(pos+1) == 'T' || seq.getX(pos+1) == 'C');
	}

	private boolean donorConstraintMinus(InputSequence<? extends Character> seq, int pos) {
		return (pos >= 2) && (seq.getX(pos-2) == 'A' || seq.getX(pos-2) == 'G') && seq.getX(pos-1) == 'C';
	}

	private boolean acceptorConstraintPlus(InputSequence<? extends Character> seq, int pos) {
		return (pos > 1) && seq.getX(pos-2) == 'A' && seq.getX(pos-1) == 'G';
	}

	private boolean acceptorConstraintMinus(InputSequence<? extends Character> seq, int pos) {
		return (seq.length() > pos + 1) && seq.getX(pos) == 'C' && seq.getX(pos+1) == 'T';
	}

	//////////////////////////////////
	
	
	private boolean stopEdgeConstraintPlus(InputSequence<? extends Character> seq, int pos) {
		if(pos < (seq.length()-2) && seq.getX(pos) == 'T') {
			return (seq.getX(pos+1) == 'A' && (seq.getX(pos+2) == 'G' || seq.getX(pos+2) == 'A'))
					|| (seq.getX(pos+1) == 'G' && seq.getX(pos+2) == 'A');
		}
		return false;
	}

	private boolean stopEdgeConstraintMinus(InputSequence<? extends Character> seq, int pos) {
		if(pos>=3 && seq.getX(pos-1) == 'A') {
			return (seq.getX(pos-2) == 'T' && (seq.getX(pos-3) == 'C' || seq.getX(pos-3) == 'T'))
					|| (seq.getX(pos-2) == 'C' && seq.getX(pos-3) == 'T');
		}
		return false;
	}
	
	/////////////////////////////////////////
	
	private boolean stopNodeConstraintPlus(InputSequence<? extends Character> seq, int pos) {
		if(pos >= 2 && seq.getX(pos-2) == 'T') {
			return (seq.getX(pos-1) == 'A' && (seq.getX(pos) == 'G' || seq.getX(pos) == 'A'))
					|| (seq.getX(pos-1) == 'G' && seq.getX(pos) == 'A');
		}
		return false;
	}

	private boolean stopNodeConstraintMinus(InputSequence<? extends Character> seq, int pos) {
		if(pos<(seq.length()-2) && seq.getX(pos+2) == 'A') {
			return (seq.getX(pos+1) == 'T' && (seq.getX(pos) == 'C' || seq.getX(pos) == 'T'))
					|| (seq.getX(pos+1) == 'C' && seq.getX(pos) == 'T');
		}
		return false;
	}
	
	@Override
	public CacheStrategySpec getCacheStrategy() {
		return new CacheStrategySpec(CacheStrategy.UNSPECIFIED);
	}

}
