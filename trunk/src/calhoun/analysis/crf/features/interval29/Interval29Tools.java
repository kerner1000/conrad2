package calhoun.analysis.crf.features.interval29;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.util.Assert;

public class Interval29Tools {
	private static final Log log = LogFactory.getLog(Interval29Tools.class);

	static protected enum Constraint {NONE, NEVER, PSTART, PDON, PACC, PSTOP, MSTART, MDON, MACC, MSTOP, PCODE, MCODE, PKEEPE, PKEEPI, MKEEPE, MKEEPI, PSTOPPED, MSTARTED, PWILLSTART, MWILLSTOP};
	static protected Constraint[] edgeConstraints;
	static protected Constraint[] nodeConstraints;
	static int numStates;

	static {
		log.debug("Setting up constraints in Interval29Tools");
		
		numStates = 29;
		
		// Setup the node constraints
		nodeConstraints = new Constraint[numStates];
		for (int j=0; j<numStates; j++) {
			nodeConstraints[j] = Constraint.NONE;
		}
		nodeConstraints[1] = Constraint.PCODE;
		nodeConstraints[2] = Constraint.PCODE;
		nodeConstraints[3] = Constraint.PCODE;
		nodeConstraints[7] = Constraint.MCODE;
		nodeConstraints[8] = Constraint.MCODE;
		nodeConstraints[9] = Constraint.MCODE;
		log.debug("The node constraints are as follows:");
		for (int i=0; i<numStates; i++) {
			log.debug("  " + i + "  --  " + nodeConstraints[i]);
		}
		
		// setup the edge constraints
		// edgeConstraint[i] is the constraint on the edge
		//   from state floor(i/numStates)
		//   to state i%numStates
		edgeConstraints = new Constraint[numStates*numStates];
		
		// The transition is impossible except when explicitly allowed below
		for(int i=0; i<numStates; ++i) {
			for(int j=0; j<numStates; ++j) {
				edgeConstraints[i*numStates + j] = Constraint.NEVER;
			}
		}
		
		// By default, self-transitions are allowed
		for(int i=0; i<numStates; i++) {
			edgeConstraints[i*numStates + i] = Constraint.NONE;			
		}
		// e-ig -> integenic
		edgeConstraints[(14)*numStates + (0)] = Constraint.PSTOPPED;
		// em-ig -> intergenic
		edgeConstraints[(22)*numStates + (0)] = Constraint.MSTARTED;
		for(int i=0; i<3; i++) {
			// intergenic -> ig-e
			edgeConstraints[(0)*numStates + (i+13)] = Constraint.PWILLSTART; //Constraint.PSTART;
			// intergenic -> ig-em
			edgeConstraints[(0)*numStates + (i+21)] = Constraint.MWILLSTOP; //Constraint.MSTOP;
			
			// Put constraints on EXON SIDE of intergenic-exon boundaries
			// ig-e_-> e_i
			edgeConstraints[(13)*numStates + (i+1)] = Constraint.PSTART;
			// ig-em -> e_im
			edgeConstraints[(21)*numStates + (i+7)] = Constraint.MSTOP;			
			// e_i -> e-ig
			edgeConstraints[(i+1)*numStates + (14)] = Constraint.PSTOP;
			// e_im -> em-ig
			edgeConstraints[(i+7)*numStates + (22)] = Constraint.MSTART;
			
			// exon-exon
			// e_i -> e_i
			edgeConstraints[(i+1)*numStates + (i+1)] = Constraint.PCODE;
			// e_im -> e_im
			edgeConstraints[(i+7)*numStates + (i+7)] = Constraint.MCODE;
					
			// Put constraints on BOTH SIDES of intron-exon boundaries
			for(int j=0; j<3; j++) {
				// e_i -> e-i_j
				edgeConstraints[(i+1)*numStates + (j+15)] = Constraint.PDON;
				// i_i -> i-e_j
				edgeConstraints[(i+4)*numStates + (j+18)] = Constraint.PACC;
				// e_im -> em-i_jm
				edgeConstraints[(i+7)*numStates + (j+23)] = Constraint.MACC;
				// i_im -> im-e_jm
				edgeConstraints[(i+10)*numStates + (j+26)] = Constraint.MDON;
			}
			// e-i_i -> i_i (intron_i, abbr.)
			edgeConstraints[(i+15)*numStates + (i+4)] = Constraint.PKEEPI;
			// i-e_i -> e_i
			edgeConstraints[(i+18)*numStates + (i+1)] = Constraint.PKEEPE;
			// em-i_im -> i_im 
			edgeConstraints[(i+23)*numStates + (i+10)] = Constraint.MKEEPI;
			// im-e_im -> e_im
			edgeConstraints[(i+26)*numStates + (i+7)] = Constraint.MKEEPE;
		}
		
//		log.warn("The transition constraints are as follows:");
//		for (int i=0; i<numStates; i++) {
//			String s = "";
//			for (int j=0; j<numStates; j++) {
//				s += edgeConstraints[i*numStates + j] + "\t";
//			}
//			System.out.println(s);
//			System.out.println("");
//		}
	}
	
	static protected int check012(int x) {
		Assert.a(x>=0, "x is " + x);
		Assert.a(x<=2, "x is " + x);
		return x;
	}
	
	static protected void verify(ModelManager modelInfo) {
		Assert.a(modelInfo.getNumStates()==29);
		
		Assert.a(modelInfo.getStateName(0).equals("intergenic"));
		Assert.a(modelInfo.getStateName(1).equals("exon0"));
		Assert.a(modelInfo.getStateName(2).equals("exon1"));
		Assert.a(modelInfo.getStateName(3).equals("exon2"));
		Assert.a(modelInfo.getStateName(4).equals("intron0"));
		Assert.a(modelInfo.getStateName(5).equals("intron1"));
		Assert.a(modelInfo.getStateName(6).equals("intron2"));
		Assert.a(modelInfo.getStateName(7).equals("exon0m"));
		Assert.a(modelInfo.getStateName(8).equals("exon1m"));
		Assert.a(modelInfo.getStateName(9).equals("exon2m"));
		Assert.a(modelInfo.getStateName(10).equals("intron0m"));
		Assert.a(modelInfo.getStateName(11).equals("intron1m"));
		Assert.a(modelInfo.getStateName(12).equals("intron2m"));
		// XXX: add Asserts for rest of states
	}
	
	static protected List<TrainingSequence<?>> checkValidTransitions(List<? extends TrainingSequence<?>> data) {
		List<TrainingSequence<?>> goodData = new ArrayList<TrainingSequence<?>>();
		for(TrainingSequence<?> seq : data) {
			boolean validSequence = true;
			for (int pos=1; pos<seq.length(); pos++) { // note start at one not zero, so can look back at prevState
				int state = seq.getY(pos);
				int prevState = seq.getY(pos-1);
				if (Interval29Tools.edgeConstraints[prevState*Interval29Tools.numStates + state] == Interval29Tools.Constraint.NEVER) {
					System.out.println("bad: " + prevState + " " + state);
					validSequence = false;
					//Assert.a(false,"pos = "+pos+" prevState = " + modelInfo.getStateName(prevState) + "   State = " + modelInfo.getStateName(state));  // A nice side effect of making sure the input sequence is legal, can omit this if you want to.
					break;
				}
			}
			if (validSequence) {
				goodData.add(seq);
			}
		}
		return goodData;
	}
}
