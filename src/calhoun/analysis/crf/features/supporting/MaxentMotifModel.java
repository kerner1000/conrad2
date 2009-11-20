package calhoun.analysis.crf.features.supporting;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.statistics.BasicStats;
import calhoun.util.Assert;

public class MaxentMotifModel  {
	private static final long serialVersionUID = -7659288739348604129L;
	private static final Log log = LogFactory.getLog(MaxentMotifModel.class);
	boolean debug = log.isDebugEnabled();
	
	// This class is written to have a single publicly static available function
	// "trainMaxentDistributionUsingAllPairwiseConstraints", and everything else
	// is private.
	//
	// This motifs are assumed to all be of same length and at each position take
	// a value between 0-3.
	
	public static double[] trainMaxentDistributionUsingAllPairwiseConstraints(List<int[]> motifExamples, int span, int nIter, double pseudocount) {
		
		int nExamples = motifExamples.size();
		if (nExamples == 0) { log.warn("Warning -- attempting to train a maxent distribution without any examples; the flat distribution will eventually be returned."); }
		for (int j=0; j<nExamples; j++) {
			int[] motif = motifExamples.get(j);
			Assert.a(motif.length == span);
			for (int k=0; k<span; k++) {
				Assert.a( (motif[k]>=0) && (motif[k]<4) );
			}
		}
		
		List<Constraint> motifConstraints = makeAllPairwiseConstraints(motifExamples,span, pseudocount);
		
		log.debug("The numebr of motifConstraints is " + motifConstraints.size());
		
		double[] ret = trainMaxentDistribution(motifConstraints,span, nIter);
		
		return ret;
	}
	
	
	private static List<Constraint> makeAllPairwiseConstraints(List<int[]> motifExamples, int motifLen, double pseudocount) {
		int nMotif = motifExamples.size();
		Assert.a(motifLen > 1);
		for (int i=1; i<nMotif; i++) { Assert.a(motifExamples.get(i).length == motifLen); }
		
		List<Constraint> ret = new ArrayList<Constraint>();
		
		for (int pos1 = 0; pos1<(motifLen-1); pos1++) {
			for (int pos2=pos1+1; pos2<motifLen; pos2++) {
				double[] counts = new double[16];
				for (int i=0; i<16; i++) { counts[i] = 0.0; }
				
				for (int j=0; j<nMotif; j++) {
					int hash = 4*motifExamples.get(j)[pos1] + motifExamples.get(j)[pos2];
					counts[hash] += 1.0;
				}
				
				double total = 0.0;
				for (int i=0; i<16; i++) { total += counts[i]; }
				for (int i=0; i<16; i++) { counts[i] /= total; }
			
				Constraint c = new Constraint(motifLen, pos1,pos2, motifExamples, pseudocount);
				ret.add( c );
				
			}
		}
		
		return ret;
	}

	private static double[] trainMaxentDistribution( List<Constraint> motifConstraints  ,  int span, int nIter ) {
		int 	hSize = 1; for (int j=0; j<span; j++) { hSize *= 4; }
		
		double[] ret = new double[hSize];
		for (int j=0; j<hSize; j++) {
			ret[j] = 1/((double) hSize);
		}
		
		int nCon = motifConstraints.size();
		if ( nCon == 0 ) {
			log.warn("Warning -- no constraints, returning maximum entropy distribution");
			return ret;
		}
		
		for (int iter=0; iter<nIter; iter++) {
			int cNum = (int) (nCon*Math.random());
			log.debug("Enforcing constrain number " + cNum + " which is " + motifConstraints.get(cNum).stringSummary() );	
			ret = motifConstraints.get(cNum).enforce(ret);
		}
		
		return ret;		
	}
	
	private static class Constraint {
		int span=-1;
		int pos1;
		int pos2;
		double[] prob;

		int msize,size;
		static int[] sixteen, newsixteen, many, newmany;
		
		public Constraint(int newspan, int pos1, int pos2, List<int[]> motifExamples, double pseudocount ) {

			if (span != newspan) {
				span = newspan;
				msize=1; for (int j=0; j<(span-2); j++) { msize *= 4; }
				size = 16*msize;
			
				sixteen = new int[16];
				newsixteen = new int[16];
				many = new int[msize];
				newmany = new int[msize];
			}	
			
			this.pos1 = pos1;
			this.pos2 = pos2;
					
			Assert.a(0<=pos1);
			Assert.a(pos1<pos2);
			Assert.a(pos2<span);
			Assert.a(2<=span);
			
			train(motifExamples, pseudocount);
		}

		
		private void train(List<int[]> motifExamples, double pseudocount) {
			
			int nMotif = motifExamples.size();
			for (int i=0; i<nMotif; i++) { Assert.a(motifExamples.get(i).length == span); }
			
			prob = new double[16];
			for (int i=0; i<16; i++) { prob[i] = pseudocount; }
			
			for (int j=0; j<nMotif; j++) {
				int hash = 4*motifExamples.get(j)[pos1] + motifExamples.get(j)[pos2];
				prob[hash] += 1.0;
			}
			
			double total = 0.0;
			for (int i=0; i<16; i++) { total += prob[i]; }
			for (int i=0; i<16; i++) { prob[i] /= total; }
			
			double sum = BasicStats.sumDoubleArray(prob);
			Assert.a( (sum > 0.999) && (sum < 1.001) );
		}

		
		public double[] enforce(double[] pp) {

			
			double[] qq = pp;
			
			
			/* We wish first to represent the range from 0 to size-1 as the cross-sum of two
			 * arrays of integers: one representing the 16 probabilities constrained by the constraint,
			 * and one representing the values of the big joint that are all within the same cell
			 * of the constraint.
			 */
			
			many[0]=0;        int nM=1;
			sixteen[0]=0;     int nS=1;
			
			for (int pos=0; pos<span; pos++) {
				if ( (pos==pos1) || (pos==pos2) ) {
					//List<Integer> newsixteen = new ArrayList<Integer>(); 
					for (int j=0; j<nS; j++) {
						int temp = 4*sixteen[j];
						newsixteen[4*j]   = temp;
						newsixteen[4*j+1] = temp+1;
						newsixteen[4*j+2] = temp+2;
						newsixteen[4*j+3] = temp+3;
					}
					nS*=4;
					for (int j=0; j<nS; j++) { sixteen[j] = newsixteen[j]; }
					//List<Integer> newmany = new ArrayList<Integer>(); 
					for (int j=0; j<nM; j++) {
						many[j] *= 4;
					}
				} else {
					//List<Integer> newsixteen = new ArrayList<Integer>(); 
					for (int j=0; j<nS; j++) {
						sixteen[j] *= 4;
					}
					//List<Integer> newmany = new ArrayList<Integer>(); 
					for (int j=0; j<nM; j++) {
						int temp = 4*many[j];
						newmany[4*j] = temp;
						newmany[4*j+1] = temp+1;
						newmany[4*j+2] = temp+2;
						newmany[4*j+3] = temp+3;
					}
					nM *= 4;
					for (int j=0; j<nM; j++) { many[j] = newmany[j]; }
				}	
			}
			
			
			Assert.a(nS == 16);
			Assert.a(nM == msize);
			
			
			double changeneeded = 0.00;
			for (int z=0; z<16; z++) {
				int base = sixteen[z];
				double total = 0.0;
				for (int t=0; t<nM; t++) {
					total += qq[base + many[t]];
				}
				if ( prob[z]>0 ) {
					Assert.a(total>0);
					changeneeded += Math.abs(prob[z] - total);
					double mult = prob[z]/total;
					for (int t=0; t<nM; t++) {
						int temp = base + many[t]; 
						qq[temp] = qq[temp]*mult;
					}
				}
			}

			System.out.println("amount of change needed to enforce constraint was " + changeneeded );
			
			return qq;
		}

		public String stringSummary() {
			String ret = "constraint_pos1=" + pos1 + "_pos2="+pos2 + "_span="+span;
			return ret;
		}
	}
	
}

