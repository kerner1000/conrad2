package calhoun.analysis.crf.features.interval13;

import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

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
public class GeneConstraintsInterval13  extends AbstractFeatureManager<Character> implements FeatureManagerEdge<Character>, FeatureManagerNode<Character> {
	@SuppressWarnings("unused")
	private static final Log log = LogFactory.getLog(GeneConstraintsInterval13.class);
	private static final long serialVersionUID = 3041359216265032511L;
	
	public String getFeatureName(int featureIndex) {
		return "Gene constraints for the model Interval13";
	}

	/** This is a constraint class, so we don't return features */
	public int getNumFeatures() {
		return 0;
	}

	public void train(int startingIndex, ModelManager modelInfo, List<? extends TrainingSequence<? extends Character>> data) {
		Interval13Tools.verify(modelInfo);
	}
	
	public void evaluateEdge(InputSequence<? extends Character> seq, int pos, int prevState, int state, FeatureList result) {
		boolean valid = true;
		
		int eind,iind;
		
		switch(Interval13Tools.edgeConstraints[prevState*Interval13Tools.numStates + state]) {
		case NONE:
			break;
		case NEVER:
			Assert.a(false);
			break;
		case PSTART:
			eind = Interval13Tools.check012(state-1);
			if ((pos-eind)%3 != 0) { valid = false; break; } 
			valid = startConstraintPlus(seq, pos);
			break;
		case PDON:
			iind = Interval13Tools.check012(state-4);
			eind = Interval13Tools.check012(prevState-1);
			if ((pos-eind+iind)%3 != 0) { valid = false; break; } 
			valid = donorConstraintPlus(seq, pos);
			break;
		case PACC:
			iind = Interval13Tools.check012(prevState-4);
			eind = Interval13Tools.check012(state-1);
			if ((pos-eind+iind)%3 != 0) { valid = false; break; } 
			valid = acceptorConstraintPlus(seq, pos);
			break;
		case PSTOP:
			eind = Interval13Tools.check012(prevState-1);
			if ((pos-eind)%3 != 0) { valid = false; break; } 
			valid = stopEdgeConstraintPlus(seq, pos);
			break;
		case MSTART:
			eind = Interval13Tools.check012(prevState-7);
			if ((pos-eind)%3 != 0) { valid = false; break; } 
			valid = startConstraintMinus(seq, pos);
			break;
		case MDON:
			iind = Interval13Tools.check012(prevState-10);
			eind = Interval13Tools.check012(state-7);
			if ((pos-eind-iind)%3 != 0) { valid = false; break; } 
			valid = donorConstraintMinus(seq, pos);
			break;
		case MACC:
			iind = Interval13Tools.check012(state-10);
			eind = Interval13Tools.check012(prevState-7);
			if ((pos-eind-iind)%3 != 0) { valid = false; break; } 			
			valid = acceptorConstraintMinus(seq, pos);
			break;
		case MSTOP:
			eind = Interval13Tools.check012(state-7);
			if ((pos-eind)%3 != 0) { valid = false; break; } 
			valid = stopEdgeConstraintMinus(seq, pos);
			break;
		case PCODE: // redundant with node invalidation below
			eind = Interval13Tools.check012(state-1);
			if ( (pos-eind)%3 == 2) {
				valid = !stopNodeConstraintPlus(seq, pos);
			}
			break;
		case MCODE: // redundant iwth node evaluation below
			eind = Interval13Tools.check012(state-7);
			if ( (pos-eind)%3==0) {
				valid = !stopNodeConstraintMinus(seq, pos);
			}
			break;
		default:
			Assert.a(false);
		}
		
		// This debugging code is pretty clutch
//		if(!valid) {
//			String str = "";
//			for (int i = -8; i < 9; i++) {
//				str += seq.getX(pos+i);
//			}
//			System.out.println("        v        ");
//			System.out.println(str);
//			System.out.println(seq.toString());
//		}			
		
		if(valid == false)
			result.invalidate();
	}

	public void evaluateNode(InputSequence<? extends Character> seq, int pos, int state, FeatureList result) {
		boolean valid = true;
	
		int eind;
		
		switch(Interval13Tools.nodeConstraints[state]) {
		case NONE:
			break;
		case NEVER:
			Assert.a(false);
			break;
		case PCODE:
			eind = Interval13Tools.check012(state-1);
			if ( (pos-eind)%3 == 2) {
				valid = !stopNodeConstraintPlus(seq, pos);
			}
			break;
		case MCODE:
			eind = Interval13Tools.check012(state-7);
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
		boolean ret = (pos >= 2) && (seq.getX(pos-2) == 'A' || seq.getX(pos-2) == 'G') && seq.getX(pos-1) == 'C';
		//if(!ret) log.warn("Seq wrong at MDON");
		return ret;
	}

	private boolean acceptorConstraintPlus(InputSequence<? extends Character> seq, int pos) {
		return (pos > 1) && seq.getX(pos-2) == 'A' && seq.getX(pos-1) == 'G';
	}

	private boolean acceptorConstraintMinus(InputSequence<? extends Character> seq, int pos) {
		boolean ret = (seq.length() > pos + 1) && seq.getX(pos) == 'C' && seq.getX(pos+1) == 'T';
		//if(!ret) log.warn("Seq wrong at MACC - expected CT but was "+seq.getX(pos)+seq.getX(pos+1));
		return ret;
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
			boolean ret = (seq.getX(pos-2) == 'T' && (seq.getX(pos-3) == 'C' || seq.getX(pos-3) == 'T'))
					|| (seq.getX(pos-2) == 'C' && seq.getX(pos-3) == 'T');
			//if(!ret) log.warn("Seq wrong at MSTOP edge entry");
			return ret;
		}
		//log.warn("Seq wrong at MSTOP edge exit");
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
			boolean ret = (seq.getX(pos+1) == 'T' && (seq.getX(pos) == 'C' || seq.getX(pos) == 'T'))
					|| (seq.getX(pos+1) == 'C' && seq.getX(pos) == 'T');
			return ret;
		}
		return false;
	}
	
	@Override
	public CacheStrategySpec getCacheStrategy() {
		return new CacheStrategySpec(CacheStrategy.UNSPECIFIED);
	}

}
