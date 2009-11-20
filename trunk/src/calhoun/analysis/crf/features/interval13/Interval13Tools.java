package calhoun.analysis.crf.features.interval13;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.ModelManager;
import calhoun.util.Assert;

public class Interval13Tools {
	private static final Log log = LogFactory.getLog(Interval13Tools.class);

	static protected enum Constraint {NONE, NEVER, PSTART, PDON, PACC, PSTOP, MSTART, MDON, MACC, MSTOP, PCODE, MCODE};
	static protected Constraint[] edgeConstraints;
	static protected Constraint[] nodeConstraints;
	static int numStates;

	static {
		log.debug("Setting up constraints in Interval13Tools");
		
		numStates = 13;
		
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
		for (int i=0; i<13; i++) {
			log.debug("  " + i + "  --  " + nodeConstraints[i]);
		}
		
		
		// setup the edge constraints
		edgeConstraints = new Constraint[numStates*numStates];
		
		// The transition is imossible except when explicitly allowed below
		for(int i=0; i<numStates; ++i) {
			for(int j=0; j<numStates; ++j) {
				edgeConstraints[i*numStates + j] = Constraint.NEVER;
			}
		}
		
		// By default, self-transitions are allowed
		for(int i=0; i<numStates; i++) {
			edgeConstraints[i*numStates + i] = Constraint.NONE;			
		}
		
		for (int i=0; i<3; i++) {
			// Following six constraints enforce the open reading frame
			edgeConstraints[(0)*numStates + (i+1)] = Constraint.PSTART;
			edgeConstraints[(i+1)*numStates + (0)] = Constraint.PSTOP;
			edgeConstraints[(i+1)*numStates + (i+1)] = Constraint.PCODE;

			edgeConstraints[(0)*numStates + (i+7)] = Constraint.MSTOP;
			edgeConstraints[(i+7)*numStates + (0)] = Constraint.MSTART;
			edgeConstraints[(i+7)*numStates + (i+7)] = Constraint.MCODE;
			
			for (int j=0; j<3; j++) {
				// Following four constraints enforce the splice rules
				edgeConstraints[(i+1)*numStates + (j+4)] = Constraint.PDON;
				edgeConstraints[(i+4)*numStates + (j+1)] = Constraint.PACC;				
				edgeConstraints[(i+10)*numStates + (j+7)] = Constraint.MDON;
				edgeConstraints[(i+7)*numStates + (j+10)] = Constraint.MACC;
			}
		}
		
		log.debug("The transition constraints are as follows:");
		for (int i=0; i<13; i++) {
			String s = "";
			for (int j=0; j<13; j++) {
				s += edgeConstraints[i*numStates + j] + "\t";
			}
			log.debug(s);
		}
	}
	
	static protected int check012(int x) {
		Assert.a(x>=0);
		Assert.a(x<=2);
		return x;
	}
	
	static protected void verify(ModelManager modelInfo) {
		Assert.a(modelInfo.getNumStates()==13);
		
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
	}
}
