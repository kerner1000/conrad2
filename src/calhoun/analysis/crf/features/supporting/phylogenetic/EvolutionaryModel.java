package calhoun.analysis.crf.features.supporting.phylogenetic;

import java.io.Serializable;
import java.util.Arrays;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.io.MultipleAlignmentInputSequence.MultipleAlignmentColumn;
import calhoun.seq.KmerHasher;
import calhoun.util.Assert;
import cern.colt.matrix.DoubleMatrix2D;

public class EvolutionaryModel implements Serializable {
	private static final long serialVersionUID = -4531626481183209673L;

	private static final Log log = LogFactory.getLog(EvolutionaryModel.class);

	// Fundamental information:
	PhylogeneticTreeFelsensteinOrder             T;
	double[]                     pi;
	NucleotideSubstitutionModel  R;
	int numSpecies;

	
	// Derived information, or precomputed/reserved for efficiency:
	int[] ileft,iright;  // indices of left and right child nodes
	double[][] Tleft;  //Transition matrices for the branches going left
	double[][] Tright;  // transition matrices for branches going right
	double[][] P;   // space in which Felsenstein algorithm recursions will be performed.

	static KmerHasher hforward = new KmerHasher(KmerHasher.ACGTother,1);
	static KmerHasher hbackward = new KmerHasher(KmerHasher.ACGTotherRC,1);	
	
	public EvolutionaryModel(PhylogeneticTreeFelsensteinOrder T,
			double[] pi, NucleotideSubstitutionModel  R ) {
		this.T = T;
		this.pi = pi;
		this.R = R;
		numSpecies = T.numSpecies();
		setup();
	}

	private void setup() {
		ileft = T.getileft();
		iright = T.getiright();
		double[] bleft = T.getbleft();
		double[] bright = T.getbright();
		
		Tleft  = new double[T.nSteps][];
		Tright = new double[T.nSteps][];
		for (int j=0; j<T.nSteps; j++) {
			Tleft[j] = createArrayFromTransitionMatrix(R.transitionMatrix(bleft[j]));
			Tright[j] = createArrayFromTransitionMatrix(R.transitionMatrix(bright[j]));
		}
		Assert.a(ileft.length == T.nSteps);
		Assert.a(iright.length == T.nSteps);
		
		P = new double[T.nNodes][4];
	}
	
	
	public double logprobRC(MultipleAlignmentColumn col, boolean conditionref) {
		return logprob(col,conditionref,hbackward);
	}

	public double logprob(MultipleAlignmentColumn col, boolean conditionref) {
		return logprob(col,conditionref,hforward);
	}
	
	
	private double logprob(MultipleAlignmentColumn C,boolean conditionref, KmerHasher h) {
			
		if ( C.numSpecies() != numSpecies ) {
			Assert.a(false,"C.numspecies is " + C.numSpecies() + "  and numSpecies is " + numSpecies);
		}
		
		for (int i=0; i<numSpecies; i++) {
			int x = h.hash(C.nucleotide(i));
			for (int j=0; j<4; j++) { 
				P[i][j] = x >= 4 || x == j ? 1.0 : 0.0; 
			}
		}
		
		for (int step=0; step<T.nSteps; step++) {
			int node = step + numSpecies;			
			felsenstein(P[ileft[step]],Tleft[step],P[iright[step]],Tright[step],P[node]);
		}
		double prob = 0;
		for (int i=0; i<4; i++) {
			prob += pi[i] * P[T.numNodes()-1][i];
		}
		
		if (conditionref) {
			for (int i=1; i<T.numSpecies(); i++) {
				Arrays.fill(P[i], 1.0);
			}
			for (int step=0; step<T.nSteps; step++) {
				int node = step + numSpecies;			
				felsenstein(P[ileft[step]],Tleft[step],P[iright[step]],Tright[step],P[node]);
			}
			double denom = 0;
			for (int i=0; i<4; i++) {
				denom += pi[i] * P[T.numNodes()-1][i];
			}
			if ( !(prob/denom < 1.00000001) ) {
				Assert.a(false , "prob=" + prob + "  denom="+denom);
			}
			
			prob = prob/denom;
		}
		
		if (!(prob > 0)) {
			Assert.a(false,"prob="+prob);
		}
		if ( !(prob < 1.00000001) ) {
			Assert.a(false , "prob=" + prob );
		}
		return Math.log(prob);	
	}

	private static double[] createArrayFromTransitionMatrix(DoubleMatrix2D R) {
		double[] ret = new double[16];
		for(int i = 0; i<4; ++i) {
			for(int j = 0; j<4; ++j) {
				ret[i*4+j] = R.getQuick(i,j);
			}			
		}
		return ret;
	}
	
	private static void felsenstein(double[] lp, double[] lT,double[] rp, double[] rT,double[] pp) {
		for (int i=0; i<4; i++) { 
			double leftprob=0.0,  rightprob=0.0;
			for (int j=0; j<4; j++) {
				leftprob += lT[i*4 + j]*lp[j];
				rightprob += rT[i*4 + j]*rp[j];
			}
			
			pp[i] = leftprob*rightprob;
		}
		return;
	}

	public void summarize() {
		log.debug("Evolutionary model, initial probabilities:   pi = " + pi[0] + "\t" + pi[1] + "\t" + pi[2] + "\t" + pi[3]);
		R.summarize();
	}

}
