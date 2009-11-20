package calhoun.analysis.crf.features.supporting.phylogenetic;

import java.io.Serializable;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.util.Assert;
import calhoun.util.ColtUtil;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

public abstract class NucleotideSubstitutionModel implements Serializable{
	private static final Log log = LogFactory.getLog(NucleotideSubstitutionModel.class);

	/* References:
	 * 
	 * http://umber.sbs.man.ac.uk/resources/phase/manual/node65.html
	 * (a description of nucleotide substitution models implemented by the program PHASE by Gowri-Shankar Vivek
	 * 
	 * Whelan, S., P. Lio, and N. Goldman
	 * 2001. Molecular phylogenetics: state-of-the art methods for looking into the past. 
	 * TRENDS in Genetics, 17(5):262-272. 
	 */
	
	DoubleMatrix2D R = new DenseDoubleMatrix2D(4,4); // nucleotide substitution rate matrix
	
	public abstract String getEvolutionaryModelName();
	
	public abstract double[] getParameters();
	
	public abstract void setParameters(double[] parms); // responsible for setting the transition rate matrix.
	
	public DoubleMatrix2D getRateMatrix() {
		return R;
	}
	
	public  DoubleMatrix2D transitionMatrix(double t) {
		//return exp(t*R) is what we want to do, but isit really that simple?
	
		DoubleMatrix2D S = new DenseDoubleMatrix2D(4,4);
		S.assign(R);
		ColtUtil.scalarmultiply(S,t);
	
		DoubleMatrix2D T = new DenseDoubleMatrix2D(4,4);
		//ColtUtil.exponentiate_symmetric_matrix(S,T);
		ColtUtil.exponentiate_real_matrix(S,T,100);
		
		// Now let us verify that the row sums of T are all zero.
		for (int i=0; i<4; i++) {
			double sumi = 0;
			for (int j=0; j<4; j++) {
				sumi += T.getQuick(i,j);
			}
			if ((sumi<0.999) || (sumi>1.001) ) {
				log.warn("Assertion failure imminent.  The matrix S is:");
				log.warn(ColtUtil.format(S));
				log.warn("The matrix T (which should be exponential of S) is:");
				log.warn(ColtUtil.format(T));				
				Assert.a(false);
			}
		}
		
		return T;
	}

	public abstract void summarize();
	
}
