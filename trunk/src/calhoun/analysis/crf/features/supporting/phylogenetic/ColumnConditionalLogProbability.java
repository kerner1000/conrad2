package calhoun.analysis.crf.features.supporting.phylogenetic;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.analysis.crf.io.MultipleAlignmentInputSequence.MultipleAlignmentColumn;
import calhoun.seq.KmerHasher;
import calhoun.util.Assert;
import flanagan.math.Minimisation;
import flanagan.math.MinimisationFunction;

public class ColumnConditionalLogProbability implements Serializable {
	private static final long serialVersionUID = 5837909206967310115L;

	private static final Log log = LogFactory.getLog(ColumnConditionalLogProbability.class);

	InputSequence<MultipleAlignmentColumn> lastSeq; 
	int lastPos;
	double[] vals;	
	
	List<int[]> clusters;             // will be an input
	int[] state2cluster;              // will be derived
	List<EvolutionaryModel> emodels;  // this gets trained
	private ModelManager model;
	static KmerHasher h = new KmerHasher(KmerHasher.ACGTother,1);
	Map<Integer, double[]> columnCache = new HashMap<Integer, double[]>();
	private int eModelNum;  // 0 = Kimura80Model
                            // 1 = HKY85Model
		
	public ColumnConditionalLogProbability(List<int[]> clusters, int eModelNum) {
		this.clusters = clusters;
		this.eModelNum = eModelNum;
	}
	
	public int numClusters() {
		return clusters.size();
	}

	public int state2cluster(int state) {
		return state2cluster[state];
	}

	public double condLogProb(InputSequence<? extends MultipleAlignmentColumn> seq, int pos, int state) {
		Assert.a(state < model.getNumStates());

		MultipleAlignmentColumn col = seq.getX(pos);
		return emodels.get(state2cluster[state]).logprob(col,true);
		
		/*
		if(seq != lastSeq || pos != lastPos) {
			lastSeq = seq;
			lastPos = pos;

			MultipleAlignmentColumn col = seq.getX(pos);
			Integer hash = col.getUniqueHash();
			if(columnCache.containsKey(hash)) {
				vals = columnCache.get(hash); 
			}
			else {
				for (int i=0; i<clusters.size(); i++) {
					vals = new double[clusters.size()];
					vals[i] = emodels.get(i).logprob(col,true);
				}
				columnCache.put(hash, vals);
			}
		}
		return vals[state2cluster[state]];
		*/
	}
	
	public void train( ModelManager modelInfo, final List<? extends TrainingSequence<? extends MultipleAlignmentColumn>> data) {
		model = modelInfo;
		vals = new double[clusters.size()];
		final PhylogeneticTreeFelsensteinOrder felsOrder = data.get(0).getX(0).getMultipleAlignment().getFelsensteinOrder();

		
		// Step 1: combinatorially invert clusters to get state2cluster
		int nStates = model.getNumStates();
		state2cluster = new int[nStates];
		for (int i=0; i<nStates; i++) { state2cluster[i] = -1; }

		for (int j=0; j<clusters.size(); j++) {
			int[] C = clusters.get(j);
			for (int k=0; k<C.length; k++) {
				state2cluster[C[k]] = j;
			}
		}
		
		for (int i=0; i<nStates; i++) {
			Assert.a(state2cluster[i] >= 0);
		}	
		
		
		// Step 2: Train the evolutionary models
		emodels = new ArrayList<EvolutionaryModel>();
		for (int j=0; j<clusters.size(); j++) {
			final int cluster = j;

			for (int k=0; k<clusters.get(j).length; k++) {
				log.debug("Training evolutionary model for:  " + model.getStateName(clusters.get(j)[k]));
			}
			
			// RootedBinaryPhylogeneticTree rt = 
						
			//final PhylogeneticTreeFelsensteinOrder T = new PhylogeneticTreeFelsensteinOrder();
			
			// Estimate pi based on the nucleotide frequencies in the reference sequence
			final double[] pi = new double[]{1.0,1.0,1.0,1.0};
			for(TrainingSequence<? extends MultipleAlignmentColumn> aln : data) {
				int len = aln.length();
					
				for (int ix=0; ix<len; ix++) {
					int y = aln.getY(ix);
					
					if (state2cluster[y] == cluster) {
						int x = h.hash(aln.getX(ix).nucleotide(0));
						if (x<4) { pi[x] += 1.0; }								
					}
				}
			}
			double total = pi[0] + pi[1] + pi[2] + pi[3];
			pi[0]/=total; pi[1]/=total; pi[2]/=total; pi[3]/=total;
			
			if (eModelNum==0) { // Kimura80Model
				MinimisationFunction mFunc = new MinimisationFunction() {
					public double function(double[] d) {
						double[] ed = new double[2];
						ed[0] = Math.exp(d[0]);
						ed[1] = Math.exp(d[1]);
						
						Kimura80Model R = new Kimura80Model(ed);
						EvolutionaryModel M = new EvolutionaryModel(felsOrder,pi,R);
						
						double ret = 0;
						for(TrainingSequence<? extends MultipleAlignmentColumn> aln : data) {
							int len = aln.length();							
							for (int ix=0; ix<len; ix++) {
								int y = aln.getY(ix);
								if (state2cluster[y] == cluster) {
									ret += M.logprob(aln.getX(ix),false);									
								}
							}
						}					
						return -ret;
					}
				};						
				
				// The standard mantra for minimizing the function mFunc defined above 
				int maxIter = 50;
				final int nParm   = 2;
				Minimisation m = new Minimisation();
				m.setNmax(maxIter);
				double[] starts = new double[nParm];
				Arrays.fill(starts, 0.1);
				double[] steps = new double[nParm];
				Arrays.fill(steps, 0.1);		
				m.nelderMead(mFunc, starts, steps);
				if(!m.getConvStatus()) {
					log.warn("WARNING - Nelder-Mead routine says convergence was not reached");
					// throw new ErrorException("Convergence not reached.");
				}
				double[] results = m.getParamValues();
				double[] eresults = new double[]{Math.exp(results[0]),Math.exp(results[1])};
				
				emodels.add(new EvolutionaryModel(felsOrder,pi,new Kimura80Model(eresults)));		
				emodels.get(cluster).summarize();
			} else if (eModelNum==1) { // HKY85Model
				MinimisationFunction mFunc = new MinimisationFunction() {
					public double function(double[] d) {
						double[] ed = new double[5];
						ed[0] = Math.exp(d[0]);
						ed[1] = Math.exp(d[1]);
						ed[2] = pi[0];
						ed[3] = pi[1];
						ed[4] = pi[2];
						
						HKY85Model R = new HKY85Model(ed);
						EvolutionaryModel M = new EvolutionaryModel(felsOrder,pi,R);
						
						double ret = 0;
						for(TrainingSequence<? extends MultipleAlignmentColumn> aln : data) {
							int len = aln.length();							
							for (int ix=0; ix<len; ix++) {
								int y = aln.getY(ix);
								if (state2cluster[y] == cluster) {
									ret += M.logprob(aln.getX(ix),false);									
								}
							}
						}					
						return -ret;
					}
				};						
				
				// The standard mantra for minimizing the function mFunc defined above 
				int maxIter = 50;
				final int nParm   = 2;
				Minimisation m = new Minimisation();
				m.setNmax(maxIter);
				double[] starts = new double[nParm];
				Arrays.fill(starts, 0.1);
				double[] steps = new double[nParm];
				Arrays.fill(steps, 0.1);		
				m.nelderMead(mFunc, starts, steps);
				if(!m.getConvStatus()) {
					log.warn("WARNING - Nelder-Mead routine says convergence was not reached");
					// throw new ErrorException("Convergence not reached.");
				}
				double[] results = m.getParamValues();
				double[] eresults = new double[]{Math.exp(results[0]),Math.exp(results[1]),pi[0],pi[1],pi[2]};
				
				emodels.add(new EvolutionaryModel(felsOrder,pi,new HKY85Model(eresults)));		
				emodels.get(cluster).summarize();
			} else { Assert.a(false); }
		}
		
		Assert.a(emodels.size() == clusters.size());
		log.debug("Just trained the Felsenstein Features");
	}
	
}
