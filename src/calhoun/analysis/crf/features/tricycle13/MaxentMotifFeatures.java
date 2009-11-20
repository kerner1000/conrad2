package calhoun.analysis.crf.features.tricycle13;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.AbstractFeatureManager;
import calhoun.analysis.crf.FeatureList;
import calhoun.analysis.crf.FeatureManagerEdge;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.features.supporting.MarkovPredictorLogprob;
import calhoun.analysis.crf.features.supporting.MaxentMotifModel;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.seq.KmerHasher;
import calhoun.util.Assert;

public class MaxentMotifFeatures extends AbstractFeatureManager<Character> implements FeatureManagerEdge<Character> {
	private static final long serialVersionUID = -7659288739348604129L;
	private static final Log log = LogFactory.getLog(MaxentMotifFeatures.class);
	boolean debug = log.isDebugEnabled();
	
	/* This is intended as an upgrade from PWM models of motifs such as splice sites, to 
	 * a more refined model callex maximum entropy model (MEM) that out to have greater
	 * sensitivity and specificity for finding splice sites, and hence also lead to more
	 * accurate gene predictons.  Relative to PWMs, MEMs can capture dependencies
	 * between postions, even non-adjacent positions.
	 * 
	 * GEOMETRY:
	 * To specify a PWM, you specify its GEOMETRY, ie its span and where that span begins
	 * relative to the transition itself (the offset), and the two hidden states before and
	 * after the transition.
	 * 
	 * DOUBLECOUNTING CORRECTION:
	 * Optionally, one may wish to subtract a double-counting correction if the observed
	 * sequence in the window would have been predicted by something else.  For this you 
	 * must also specify the predictor that would have been used by default for each base
	 * (depending on the hidden state), and the sequence of hidden states to which this
	 * predictor would have been applied over the span of the PWM.
	 * At the moment, the double-counting correction is no longer optional; it is now required.
	 * 
	 * Each motif model models a particular transition; e.g. exon3->intron3 has a model
	 * different than exon2->intron2.
	 * 
	 * The MaxEnt constraint definition and iterative scaling method are being ported from
	 * Matlab, where Jade Vinson first implemented this based on paper byBurge and Yeo.
	 */
	
	
	// MEMBER VARIABLES //////////////////////////////////////////////////////
	
	// ADMINISTRATIVE OVERHEAD
	int startIx;  // The index of the first feature managed by this FeatureManager
	ModelManager model;
	KmerHasher   h; // for a single letter
	//DenseBooleanMatrix2D[] transitions;
	boolean tieFlag = false;		
	
	// GEOMETRY
	int nFeatures;	// derived from geometry
	int spanLimit = 12;
	int[] span;		// derived from geometry
	int[] offset;	// derived from geometry
	int[] tranfrom;
	int[] tranto;
	
	// DOUBLECOUNTING CORRECTION
	boolean dcflag;
	List<int[]> dcc;  // DoubleCounting Correction
	MarkovPredictorLogprob predictorlp;
	
	// OPTIMIZATION OVERHEAD
	InputSequence<? extends Character> lastSeq;
	int lastPos;
	float[] vals;
	
	// DATA THAT GETS TRAINED
	List<double[]> listprob;
	
	
	
	//	CONSTRUCTORS AND SUPPORT //////////////////////////////////////////////////////
	
	public MaxentMotifFeatures(List<int[]> geometry, List<int[]> dccorrection, List<int[]> markovhistory) {
		setThingsUp(geometry,dccorrection,markovhistory);
	}
	
	public MaxentMotifFeatures(List<int[]> geometry, List<int[]> dccorrection, List<int[]> markovhistory, List<int[]> flags) {
		tieFlag = true;
		setThingsUp(geometry,dccorrection,markovhistory);
	}
	
	private void setThingsUp(List<int[]> geometry, List<int[]> dccorrection, List<int[]> markovhistory) {
		
		// ADMINISTRATIVE
		/* geometry, one of the inputs, has the following interpretation:
		 * For each i, geometry[i] describes the geometry of one PWM feature.
		 * 0) geometry[i][0] is the span of the PWM
		 * 1) geometry[i][1] is the offset j of the PWM, so that the feature for
		 *     position i relates the following observable and hidden states:
		 *     y_(i-1), y_i, x_(i-j), x_(i-j+1), ... , x_(i-j+span-1)
		 * 2) geometry[i][2] is yprev
		 * 3) geometry[i][3] is y  */
		
		
		
		nFeatures = geometry.size();
		vals = new float[nFeatures];
		h = new KmerHasher(KmerHasher.ACGTother, 1);	
		
		// DOUBLECOUNTING CORRECTION
		dcflag = true;
		this.predictorlp = new MarkovPredictorLogprob(markovhistory);		
		this.dcc = dccorrection;	
		
		// GEOMETRY
		span = new int[nFeatures];
		offset = new int[nFeatures];
		tranfrom = new int[nFeatures];
		tranto = new int[nFeatures];
		listprob = new ArrayList<double[]>();
		for (int i=0; i<nFeatures; i++) {
			span[i]     = geometry.get(i)[0];
			offset[i]   = geometry.get(i)[1];
			tranfrom[i] = geometry.get(i)[2];
			tranto[i]   = geometry.get(i)[3];
			Assert.a(span[i] <= spanLimit);
			int len=1; for (int j=0; j<span[i]; j++) { len *= 4; }
			double[] prob = new double[len];
			listprob.add(prob);
		}	
		
		//LOTS OF ASSERTIONS
		Assert.a(dccorrection.size()==nFeatures);		
		for (int i=0; i<nFeatures; i++) {
			Assert.a(  (offset[i]>=0) && (offset[i]<= span[i])  ); // So the span of the transition is WITHIN the span of the span of the PWM
			Assert.a(dccorrection.get(i).length == span[i]);
		}	
	}
	
	
	// BORING ADMINISTRATIVE FUNCTIONS //////////////////////////////////////////////////////
	
	public int getNumFeatures() {
		if (tieFlag) { return 1; }
		return nFeatures;
	}	
	
	public String getFeatureName(int featureIndex) {
		if (tieFlag) { return "tiedMaxentMotifModels"; }
		
		int raw = featureIndex - startIx;
		
		String ret = "MaxentMotifModels.span" + span[raw] + ".offset" + offset[raw] + ".fromState." + model.getStateName(tranfrom[raw]) + ".toState."+ model.getStateName(tranto[raw]);
		
		return ret;
	}
	
	// EVALUATION FUNCTION AND SUPPORT //////////////////////////////////////////////////////
	
	public void evaluateEdge(InputSequence<? extends Character> seq, int pos, int previousState, int state, FeatureList result) {
		if(pos == 0) {
			return;
		}
		
		if((seq != lastSeq) || (pos != lastPos)) {
			lastSeq = seq;
			lastPos = pos;
			updateVals(seq, pos);
		}
		
		if (tieFlag) {
			for (int j=0; j<nFeatures; j++) {
				if( (tranfrom[j]==previousState) && (tranto[j]==state)) {
					result.addFeature(startIx, vals[j]);
				}
			}		
		} else {
			for (int j=0; j<nFeatures; j++) {
				if( (tranfrom[j]==previousState) && (tranto[j]==state)) {
					result.addFeature(startIx + j, vals[j]);
				}
			}		
		}
	}
	
	
	
	
	
	void updateVals(InputSequence<? extends Character> seq, int ix) {
		for (int j=0; j<nFeatures; j++) {
			
			int spn = span[j];
			int offset1 = this.offset[j];
			
			// If there are no missing data in the relevant window, then return log( prob according to maxent / prob according to default )
			
			float val = 0;
			if ((ix>=offset1) && ((ix-offset1+spn)<=(seq.length())) ) {
				boolean completeInformation = true;
				for (int i=0; i<spn; i++) {
					int pos = ix - offset1 + i;
					char c = seq.getX(pos);
					if (h.hash(c) == 4) { completeInformation=false; }
				}				
				
				if (completeInformation) {
					int hash = 0;
					for (int i=0; i<spn; i++) {
						int pos = ix - offset1 + i;
						char c = seq.getX(pos);
						hash = hash*4 + h.hash(c);
					}
					double exval = listprob.get(j)[hash];
					if (exval>0) {
						val = (float) Math.log(exval);
					} else {
						log.info("Refusing to take log of zero, returning a large penalty instead.");
						val = -4000;
					}
						
					if (dcflag) {						
						for (int i=0; i<spn; i++) {
							val = val - predictorlp.logprob(dcc.get(j)[i],seq,ix-offset1+i);					
						}
					}
				}
			}
			vals[j] = val;
		}			
	}
	
	
	// TRAINING FUNCTION //////////////////////////////////////////////////////
	
	public void train(int startingIndex, ModelManager modelInfo, List<? extends TrainingSequence<? extends Character>> data) {
		startIx = startingIndex;
		model = modelInfo;
		
		log.debug("Training a maxent motif feature manager");
		
		for (int i=0; i<nFeatures; i++) {
			List<int[]> motifExamples = new ArrayList<int[]>();
			// Loop through all of the training data and record all fully present (ie no missing data) examples of the motif.
			// Use the hasher to get values between 0-3.
			
			log.debug("Training a maxent motif feature with span " + span[i]);
			
			for(TrainingSequence<? extends Character> seq : data) {
				int len = seq.length();
				
				for (int ix=offset[i]; ix<(len-span[i]+offset[i]); ix++) {
					if (ix<=0) continue;
					if (ix>=len) continue;
					int yprev = seq.getY(ix-1);
					int y = seq.getY(ix);
					
					if ( (yprev == tranfrom[i]) && (y == tranto[i]) ) {
						int[] example = new int[span[i]];
						boolean completeInformation = true;
						for (int pos=0; pos<span[i]; pos++) {
							char c = seq.getX(ix - offset[i] + pos);
							example[pos] = h.hash(c);
							if (example[pos]==4) { completeInformation = false; } 
						}
						if (completeInformation) { motifExamples.add(example); }
					}
				}	
			}
						
			//List<Constraint> motifConstraints = makeAllPairwiseConstraints(motifExamples);
			
			//listprob.set( i , trainMaxentDistribution(motifConstraints,span[i]) );
			listprob.set( i , MaxentMotifModel.trainMaxentDistributionUsingAllPairwiseConstraints(motifExamples,span[i],1000,0.01) );
		}
				
		if (dcflag) {
			predictorlp.train(data);
		}
	}
}

