package calhoun.analysis.crf.features.supporting.phylogenetic;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.util.Assert;

public class HKY85Model extends NucleotideSubstitutionModel {
	private static final long serialVersionUID = -1555350959092722736L;
	private static final Log log = LogFactory.getLog(HKY85Model.class);
	// Model for nucleotide substitution proposed by Kimura in 1980.
	// Models the overall rate of nucleotide substitution and the difference
	// between rate of transitions and transversions.
	
	double ts;
	double tv;
	double piA;
	double piC;
	double piG;
	double piT;
	
	/* HKY85 Nucleotide substitution model
	 * 1. Website:         http://www.mun.ca/biology/scarr/Transitions_vs_Transversions.htm
	 * 2. Reviw article:   Pietro Lio and Nick Goldman, "Models of Molecular Evolution and Phylogeny", Genome Research, 1998.
	 * 3. Primary article: Hasegawa, Kishino, Yano, 1985, "Dating of the human-ape splitting by a molecular clock of mitochondrial DNA", J. Molecular Evolution. 22, 160-174.
	 */
	
	public HKY85Model(double[] parms) {
		setParameters(parms);
	}

	public HKY85Model() {
		// Provide some reasonable seed values
		ts  = 0.6;
		tv  = 0.2;
		piA = 0.25;
		piC = 0.25;
		piG = 0.25;
		piT = 0.25;
	}
	
	@Override
	public String getEvolutionaryModelName() {
		return "HKY85_nucleotide_evolution_model";
	}
	
	@Override
	public double[] getParameters() {
		double[] ret = new double[6];
		ret[0] = ts;
		ret[1] = tv;
		ret[2] = piA;
		ret[3] = piC;
		ret[4] = piG;
		ret[5] = piT;
		return ret;
	}


	@Override
	public void setParameters(double[] parms) {
		Assert.a(parms.length == 5);
		for (int j=0; j<5; j++) {
			Assert.a(parms[j] > 0);
		}
		
		ts  = parms[0];
		tv  = parms[1];
		piA = parms[2];
		piC = parms[3];
		piG = parms[4];
		piT = 1 - piA - piC - piG;
		Assert.a(piT>0);
		
		// Maybe here I should be using a hashing function, but let's first just write it out;
		double[][] X = new double[4][4];
		int A=0, C=1, G=2, T=3;
		// Given an i base, X[i][j] is rate of changing to base j.
		X[A][A]=0;         X[A][C]=tv*piC;   X[A][G]=ts*piG;   X[A][T]=tv*piT;
		X[C][A]=tv*piA;    X[C][C]=0;        X[C][G]=tv*piG;   X[C][T]=ts*piT; 
		X[G][A]=ts*piA;    X[G][C]=tv*piC;   X[G][G]=0;        X[G][T]=tv*piT;
		X[T][A]=tv*piA;    X[T][C]=ts*piC;   X[T][G]=tv*piG;   X[T][T]=0;
		for (int i=0; i<4; i++) {
			for (int j=0; j<4; j++) {
				if (j!=i) {
					X[i][i] -= X[i][j];
				}
			}
		}
		
		R.assign(X);
	}

	@Override
	public void summarize() {
		log.debug("HKY85 Nucleotide substitution rate matrix; ts=" + ts + "  tv=" + tv + "  piA="+piA+"  piC="+piC+"  piG="+piG+"  piT="+piT);	
	}
}
