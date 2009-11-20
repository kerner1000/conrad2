package calhoun.analysis.crf.features.interval13;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.AbstractFeatureManager;
import calhoun.analysis.crf.CacheStrategySpec;
import calhoun.analysis.crf.FeatureList;
import calhoun.analysis.crf.FeatureManagerNode;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.CacheStrategySpec.CacheStrategy;
import calhoun.analysis.crf.features.supporting.phylogenetic.EvolutionaryModel;
import calhoun.analysis.crf.features.supporting.phylogenetic.Kimura80Model;
import calhoun.analysis.crf.features.supporting.phylogenetic.PhylogeneticTreeFelsensteinOrder;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.analysis.crf.io.MultipleAlignmentInputSequence.MultipleAlignmentColumn;
import calhoun.seq.KmerHasher;
import calhoun.util.Assert;
import flanagan.math.Minimisation;
import flanagan.math.MinimisationFunction;

public class PhylogeneticLogprobInterval13 extends AbstractFeatureManager<MultipleAlignmentColumn> implements FeatureManagerNode<MultipleAlignmentColumn> {
	private static final long serialVersionUID = -7659288739348604129L;
	private static final Log log = LogFactory.getLog(PhylogeneticLogprobInterval13.class);
	
	
	int startIx;  // The index of the first feature managed by this FeatureManager
	ModelManager model;	
	boolean multipleFeatures = false;
	
	EvolutionaryModel emodelIntergenic;      // one model for a column of aligned sequence in intergenic region
	EvolutionaryModel emodelIntronic;        // one model for intronic regions
	ArrayList<EvolutionaryModel> emodelExonic; // a model for positions 0,1,2 = (A,T,G) of a codon n a coding exon.
	
	static KmerHasher hforward = new KmerHasher(KmerHasher.ACGTother,1);    // a character hasher for forward strand
	static KmerHasher hbackward = new KmerHasher(KmerHasher.ACGTotherRC,1); // a character hasher for reverse strand
	
	///////////////////////////////////////////////////////////////////////////////
		
	public PhylogeneticLogprobInterval13() { }	  // a constructor with no arguments
	
	public int getNumFeatures() {  // there is exactly one feature
		return multipleFeatures ? 5 : 1;
	}	
	 
	public String getFeatureName(int featureIndex) {
		if(multipleFeatures) {
			String[] vals = new String[] { "Intergenic", "Exon pos.", "Intron pos.", "Exon neg.", "Intron neg."};
			int feat = featureIndex - startIx;
			String table = vals[feat];
			return table+" phylogeny";
		}
		else {
			return "PhylogeneticLogProbInterval13";
		}
	}
	
	
	public void evaluateNode(InputSequence<? extends MultipleAlignmentColumn> seq, int pos, int state, FeatureList result) {

		Assert.a(state < model.getNumStates());
		MultipleAlignmentColumn col = seq.getX(pos);
		
		double val = 0.0;
		int ephase;
		int featureOffset = Integer.MIN_VALUE;
		switch (state) {
		case 0:
			val = emodelIntergenic.logprob(col,true);
			featureOffset = 0;
			break;
		case 1:
		case 2:
		case 3:
			ephase = ((pos-state+1)%3+3)%3;  //((pos-(state-1))%3 +3)%3;
			//val = emodelExonic.get(0).logprob(col,true);
			val = emodelExonic.get(ephase).logprob(col,true);
			featureOffset = 1; // + ephase;
			break;
		case 4:
		case 5:
		case 6:
			val = emodelIntronic.logprob(col,true);
			featureOffset = 2;
			break;
		case 7:
		case 8:
		case 9:
			ephase = ((-pos+state+1)%3+3)%3;   // ((-pos+2+(state-7))%3 +3)%3;
			val = emodelExonic.get(ephase).logprobRC(col,true);
			featureOffset = 3; // + ephase;
			break;			
		case 10:
		case 11:
		case 12:
			val = emodelIntronic.logprobRC(col,true);
			featureOffset = 4;
			break;
		default:
			Assert.a(false);
		}
		
		result.addFeature(startIx + (multipleFeatures ? featureOffset : 0), val);
	}
	
	
	public void train(int startingIndex, ModelManager modelInfo, final List<? extends TrainingSequence<? extends MultipleAlignmentColumn>> data) {
		startIx = startingIndex;
		model = modelInfo;
				
		final PhylogeneticTreeFelsensteinOrder felsOrder = data.get(0).getX(0).getMultipleAlignment().getFelsensteinOrder();
		
		ArrayList<boolean[]> flagsForward  = new ArrayList<boolean[]>();
		ArrayList<boolean[]> flagsBackward = new ArrayList<boolean[]>();
		for(int seqNum=0; seqNum<data.size(); seqNum++) {
			TrainingSequence<? extends MultipleAlignmentColumn> aln = data.get(seqNum);
			int len = aln.length();
			flagsForward.add(new boolean[len]);
			flagsBackward.add(new boolean[len]);					
		}
		
		
		log.debug("Training model for intergenic regions...");
		for(int seqNum=0; seqNum<data.size(); seqNum++) {
			TrainingSequence<? extends MultipleAlignmentColumn> aln = data.get(seqNum);
			int len = aln.length();
			
			boolean[] ff = flagsForward.get(seqNum);
			Assert.a(ff.length == len);
			
			boolean[] fb = flagsBackward.get(seqNum);			
			Assert.a(fb.length == len);
			
			for (int pos=0; pos<len; pos++) {
				int y = aln.getY(pos);
				if (y == 0) {
					ff[pos] = true;
					fb[pos] = true;
				} else {
					ff[pos] = false;
					fb[pos] = false;
				}
			}		
		}		
		emodelIntergenic = trainEvolutionaryModel(felsOrder,data,flagsForward, flagsBackward);
		log.debug("Evolutionary model for intergenic regions:");
		emodelIntergenic.summarize();

		

		log.debug("Training model for intronic regions...");
		for(int seqNum=0; seqNum<data.size(); seqNum++) {
			TrainingSequence<? extends MultipleAlignmentColumn> aln = data.get(seqNum);
			int len = aln.length();
			
			boolean[] ff = flagsForward.get(seqNum);
			Assert.a(ff.length == len);
			
			boolean[] fb = flagsBackward.get(seqNum);			
			Assert.a(fb.length == len);
			
			for (int pos=0; pos<len; pos++) {
				int y = aln.getY(pos);
				if ( (y == 4) || (y == 5) || (y == 6) ) {
					ff[pos] = true;
				} else {
					ff[pos] = false;
				}
				if ( (y == 10) || (y == 11) || (y == 12) ) {
					fb[pos] = true;
				} else {
					fb[pos] = false;
				}
			}		
		}		
		emodelIntronic = trainEvolutionaryModel(felsOrder,data,flagsForward, flagsBackward);
		log.debug("Evolutionary model for intronic regions:");
		emodelIntronic.summarize();

		// 	  ephase = ((pos-state+1)%3+3)%3; for states 1,2,3
		//    ephase = ((-pos+state+1)%3+3)%3; for states 10,11,12
		
		emodelExonic = new ArrayList<EvolutionaryModel>();
		for (int phase =0; phase<3; phase++) {
			log.debug("Training model for exonic regions...");
			for(int seqNum=0; seqNum<data.size(); seqNum++) {
				TrainingSequence<? extends MultipleAlignmentColumn> aln = data.get(seqNum);
				int len = aln.length();
				
				boolean[] ff = flagsForward.get(seqNum);
				Assert.a(ff.length == len);
				
				boolean[] fb = flagsBackward.get(seqNum);			
				Assert.a(fb.length == len);
				
				for (int pos=0; pos<len; pos++) {
					int y = aln.getY(pos);
					int pstate = ((pos-phase)%3 +3)%3 + 1;
					int mstate = ((phase+pos-2)%3 + 3)%3 + 7;
					if ( y == pstate ) {
						ff[pos] = true;
					} else {
						ff[pos] = false;
					}
					if ( y==mstate ) {
						fb[pos] = true;
					} else {
						fb[pos] = false;
					}
				}		
			}		
			emodelExonic.add(trainEvolutionaryModel(felsOrder,data,flagsForward, flagsBackward));
			log.debug("Evolutionary model for intronic regions:");
			emodelExonic.get(phase).summarize();
		}
		
		log.debug("Just trained all evolutionary models");
	}
	
	private EvolutionaryModel trainEvolutionaryModel(final PhylogeneticTreeFelsensteinOrder felsOrder,
			final List<? extends TrainingSequence<? extends MultipleAlignmentColumn>> data,
			final ArrayList<boolean[]> flagsForward,
			final ArrayList<boolean[]> flagsBackward) {
		
		Assert.a(flagsForward.size() == data.size());
		Assert.a(flagsBackward.size() == data.size());
		
		// Estimate pi based on the nucleotide frequencies in the reference sequence
		final double[] pi = new double[]{1.0,1.0,1.0,1.0};
		for(int seqNum=0; seqNum<data.size(); seqNum++) {
			TrainingSequence<? extends MultipleAlignmentColumn> aln = data.get(seqNum);
			int len = aln.length();
			
			boolean[] ff = flagsForward.get(seqNum);
			Assert.a(ff.length == len);
			
			boolean[] fb = flagsBackward.get(seqNum);			
			Assert.a(fb.length == len);			
			
			for (int ix=0; ix<len; ix++) {				
				if (ff[ix]) {
					int x = hforward.hash(aln.getX(ix).nucleotide(0));
					if (x<4) { pi[x] += 1.0; }								
				}
				
				if (fb[ix]) {
					int x = hbackward.hash(aln.getX(ix).nucleotide(0));
					if (x<4) { pi[x] += 1.0; }								
				}				
			}
		}
		double total = pi[0] + pi[1] + pi[2] + pi[3];
		pi[0]/=total; pi[1]/=total; pi[2]/=total; pi[3]/=total;

		
		MinimisationFunction mFunc = new MinimisationFunction() {
			public double function(double[] d) {
				double[] ed = new double[2];
				ed[0] = Math.exp(d[0]);
				ed[1] = Math.exp(d[1]);
				
				Kimura80Model R = new Kimura80Model(ed);
				EvolutionaryModel M = new EvolutionaryModel(felsOrder,pi,R);
				
				double ret = 0;
				for(int seqNum=0; seqNum<data.size(); seqNum++) {
					TrainingSequence<? extends MultipleAlignmentColumn> aln = data.get(seqNum);
					int len = aln.length();							
					
					boolean[] ff = flagsForward.get(seqNum);
					Assert.a(ff.length == len);
					
					boolean[] fb = flagsBackward.get(seqNum);			
					Assert.a(fb.length == len);
					
					for (int ix=0; ix<len; ix++) {
						if (ff[ix]) {
							ret += M.logprob(aln.getX(ix),true);									
						}
						if (fb[ix]) {
							ret += M.logprobRC(aln.getX(ix),true);
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
		}
		double[] results = m.getParamValues();
		double[] eresults = new double[]{Math.exp(results[0]),Math.exp(results[1])};
		
		return (new EvolutionaryModel(felsOrder,pi,new Kimura80Model(eresults)) );		
	}

	@Override
	public CacheStrategySpec getCacheStrategy() {
		return new CacheStrategySpec(CacheStrategy.DENSE);
	}

	/**
	 * @return Returns the multipleFeatures.
	 */
	public boolean isMultipleFeatures() {
		return multipleFeatures;
	}

	/**
	 * @param multipleFeatures The multipleFeatures to set.
	 */
	public void setMultipleFeatures(boolean multipleFeatures) {
		this.multipleFeatures = multipleFeatures;
	}
}
