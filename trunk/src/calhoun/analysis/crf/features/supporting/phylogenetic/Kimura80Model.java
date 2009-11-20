package calhoun.analysis.crf.features.supporting.phylogenetic;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.util.Assert;

public class Kimura80Model extends NucleotideSubstitutionModel {
	private static final long serialVersionUID = 4075281397219042249L;

	private static final Log log = LogFactory.getLog(Kimura80Model.class);

	// Model for nucleotide substitution proposed by Kimura in 1980.
	// Models the overall rate of nucleotide substitution and the difference
	// between rate of transitions and transversions.
	
	double ts;
	double tv;	
	
	/* Transitions (Ts) are interchanges between pyrimidines (C  T), or between purines (A  G) 
	 * Transversions (Tv) are interchanges beween purines & pyrimidines.
	 * 
	 * There are twice as many possible transversions, however transitions occur more
	 * commonly because of the molecular mechanisms involved.
	 * 
	 * (from http://www.mun.ca/biology/scarr/Transitions_vs_Transversions.htm)
	 */
	
	public Kimura80Model(double[] parms) {
		setParameters(parms);
	}

	public Kimura80Model() {
		// Provide some reasonable seed values
		ts = 0.6;
		tv = 0.2;
	}

	@Override
	public String getEvolutionaryModelName() {
		return "Kimura80_nucleotide_evolution_model";
	}
	
	@Override
	public double[] getParameters() {
		double[] ret = new double[2];
		ret[0] = ts;
		ret[1] = tv;
		return ret;
	}


	@Override
	public void setParameters(double[] parms) {
		Assert.a(parms.length == 2);
		Assert.a(parms[0] > 0);
		Assert.a(parms[1] > 0);		
		
		ts   = parms[0];
		tv = parms[1];
				
		// Maybe here I should be using a hashing function, but let's first just write it out;
		double[][] X = new double[4][4];
		int A=0, C=1, G=2, T=3;
		X[A][A]=-2*tv-ts;  X[A][C]=tv;       X[A][G]=ts;       X[A][T]=tv;
		X[C][A]=tv;        X[C][C]=-2*tv-ts; X[C][G]=tv;       X[C][T]=ts; 
		X[G][A]=ts;        X[G][C]=tv;       X[G][G]=-2*tv-ts; X[G][T]=tv;
		X[T][A]=tv;        X[T][C]=ts;       X[T][G]=tv;       X[T][T]=-2*tv-ts;
		
		R.assign(X);
	}

	@Override
	public void summarize() {
		log.debug("Kimura80 Nucleotide substitution rate matrix; ts=" + ts + "  tv=" + tv);	
	}
}
