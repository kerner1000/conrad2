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
import calhoun.analysis.crf.features.supporting.phylogenetic.ColumnConditionalLogProbability;
import calhoun.analysis.crf.features.supporting.phylogenetic.EvolutionaryModel;
import calhoun.analysis.crf.io.CompositeInput;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.analysis.crf.io.MultipleAlignmentInputSequence.MultipleAlignmentColumn;
import calhoun.seq.KmerHasher;
import calhoun.util.Assert;
import calhoun.util.DenseBooleanMatrix2D;

public class PWM_evolution extends AbstractFeatureManager<CompositeInput> implements FeatureManagerEdge<CompositeInput> {
	private static final long serialVersionUID = -7659288739348604129L;
	private static final Log log = LogFactory.getLog(PWM_evolution.class);
	boolean debug = log.isDebugEnabled();
	
	/* PWM evolution is like the position weight matrix fetaures for modeling the
	 * boundaries between two extensive states (eg a donor site separating exons from
	 * introns).  However, this feature does more because it not only does a PWM
	 * for the reference sequence but also trains a Kimura80 model for nucleotide
	 * evolution at each position of the feature (each position gets its own model).
	
	 * Note that one must subtract double-counting corrections, since these bases would
	 * otherwise have been modeled using the exon states or the intron states.
	 * We'll use the same machinery that PWM uses for knowing what to subtract (ie what is
	 * being replaced.  However, here we must subtract out not only the prior assignment to the reference sequence (as we did with PWM)
	 * but also the prior assignemnt of conditional probability of the multiple alignment
	 * column given the reference sequence and evolutionary model.
	 * 
	 * Note similarities to both PositionWeightMatrixFeatures.java
	 *   and FelsensteinFeatures.java
	 */
	
	int startIx;  // The index of the first feature managed by this FeatureManager
	ModelManager model;
	
	// Following block are things that depend explicitly on and are calculated directly
	//  from geometry, included only for convenience. 
	int nFeatures;
	int[] span;
	int[] offset;
	int[] nTrans;


	
	DenseBooleanMatrix2D[] transitions;
	
	// New requirement: each Feture, or element of geometry, describes a single transition
	// The variable geometry is the information that is needed to initialize
	List<int[]> geometry; /* For each i, geometry[i] describes the geometry of one PWM feature.
	                       * 0) geometry[i][0] is the span of the PWM
	                       * 1) geometry[i][1] is the offset j of the PWM, so that the feature for
	                       *     position i relates the following observable and hidden states:
	                       *     y_(i-1), y_i, x_(i-j), x_(i-j+1), ... , x_(i-j+span-1)
	                       * 2) geometry[i][2] is yprev
	                       * 3) geometry[i][3] is y  */
	
	// These are the parameters that need to be trained empirically.
	List<float[][]> logprob; /* logprob[i][j] is log probability of base j at position i, i=0..(span-1), j=0..3. */

	
	
	List<int[]> dcc;  // DoubleCounting Correction
	MarkovPredictorLogprob predictorlp;
	
	// This stuff is similar to FelsensteinFeatures:
	List<int[]> clusters;             // will be an input
	List<EvolutionaryModel> emodels;  // this gets trained
	int[] state2cluster;
	static KmerHasher h = new KmerHasher(KmerHasher.ACGTother,1);
	ColumnConditionalLogProbability mo;
	boolean tieFlag;

	InputSequence<? extends CompositeInput> lastSeq;
	int lastPos;
	float[] vals;
	private int nUpdate = 0;
	

	public PWM_evolution(List<int[]> geometry, List<int[]> dccorrection, List<int[]> markovhistory, List<int[]> clusters) {
		tieFlag = false;
		PWM_evolution_setup(geometry,dccorrection, markovhistory,clusters);
	}

	public PWM_evolution(List<int[]> geometry, List<int[]> dccorrection, List<int[]> markovhistory, List<int[]> clusters, List<int[]> flags) {
		tieFlag = true;
		PWM_evolution_setup(geometry,dccorrection, markovhistory,clusters);
	}

	private void PWM_evolution_setup(List<int[]> geometry1, List<int[]> dccorrection, List<int[]> markovhistory, List<int[]> clusters1) {
		dcc = dccorrection;
		predictorlp = new MarkovPredictorLogprob(markovhistory);
		mo = new ColumnConditionalLogProbability(clusters1,0); // Zero corresponds to default Kimura80 model

		this.geometry = geometry1;
		this.clusters = clusters1;
		
		setupGeometry();	
	}

	private void setupGeometry() {
		nFeatures = geometry.size();
		span = new int[nFeatures];
		offset = new int[nFeatures];
		nTrans = new int[nFeatures];
		h = new KmerHasher(KmerHasher.ACGTN, 1);
		logprob = new ArrayList<float[][]>();
		for (int i=0; i<nFeatures; i++) {
			nTrans[i] = (geometry.get(i).length - 2)/2;
			span[i] = geometry.get(i)[0];
			offset[i] = geometry.get(i)[1];
			Assert.a(offset[i]>=0); // So the span of the transition is WITHIN the span of the span of the PWM
			float[][] lp = new float[span[i]][h.range()];
			logprob.add(lp);
		}
		
		Assert.a(geometry.size()==nFeatures);
		Assert.a(dcc.size()==nFeatures);		
		for (int j=0; j<nFeatures; j++) {
			Assert.a(nTrans[j]==1);
			Assert.a(dcc.get(j).length == span[j]);
		}
	}

	public int getNumFeatures() {
		if (tieFlag) { return 1; }
		return nFeatures;
	}	
	
	public String getFeatureName(int featureIndex) {
		int raw = featureIndex - startIx;
		int[] X = geometry.get(raw);
		
		String ret = "PWM.span" + X[0] + ".offset" + X[1];
		for (int j=2; j<X.length; j+=2) {
			ret = ret + ".(" + model.getStateName(X[j]) + "," + model.getStateName(X[j+1]) + ")"; 
		}	
		return ret;
	}


	//static int count = 0;
	public void evaluateEdge(InputSequence<? extends CompositeInput> seq, int pos, int previousState, int state, FeatureList result) {
		if(pos == 0) {
			return;
		}
		
		if(seq != lastSeq || pos != lastPos) {
			lastSeq = seq;
			lastPos = pos;
			//System.out.println("Pos: "+pos+" Seq: "+seq);
			//if(++count > 10) {
			//	throw new RuntimeException();
			//}
			updateVals(seq, pos);
		}
		
		for (int j=0; j<nFeatures; j++) {
			if(transitions[j].getQuick(previousState, state)) {
				if (tieFlag) { 
					result.addFeature(startIx, vals[j]);					
				} else {
					result.addFeature(startIx + j, vals[j]);
				}
			}
		}		
	}

	void updateVals(InputSequence<? extends CompositeInput> seq, int ix) {
		nUpdate ++;
		for (int j=0; j<nFeatures; j++) {
			// Figure out what needed for Feature j, which might have several valid y pairs, at position i
			int[] geo = geometry.get(j);
			int spn = geo[0];
			int offset1 = geo[1];
			float val = 0;
			if ((ix>=offset1) && ((ix-offset1+spn)<=(seq.length())) ) {
				InputSequence<Character>                CIS = (InputSequence<Character>) seq.getComponent("ref");
				InputSequence<MultipleAlignmentColumn>  MIS = (InputSequence<MultipleAlignmentColumn>) seq.getComponent("aln");
				
				for (int i=0; i<spn; i++) {
					int pos = ix - offset1 + i;
					char c = CIS.getX(pos);
					val = val + logprob.get(j)[i][h.hash(c)];
				}
				
				Assert.a(nTrans[j]==1);
				Assert.a(dcc.get(j).length == spn);
				
				for (int i=0; i<spn; i++) {
					// dcc.get(j)[i] is a state (a number 0-12); ix-offset+i is a position
					// This is for subtracting to correct for what was added by the feature
					// MarkovPredictorLogProb and is being replaced by the PWM
					val = val - predictorlp.logprob(dcc.get(j)[i],CIS,ix-offset1+i);
					// Below subtracts the correction for the feature FelsensteinFeatures
					// which we are now replacing with the feature PWM_evolution at this position.
					val = (float) (val - mo.condLogProb(MIS,ix-offset1+i,dcc.get(j)[i]));
				}
			}
			vals[j] = val;
		}			
	}

	public void train(int startingIndex, ModelManager modelInfo, List<? extends TrainingSequence<? extends CompositeInput>> data) {
		startIx = startingIndex;
		model = modelInfo;

		vals = new float[nFeatures];
		
		List<TrainingSequence<Character>> LTSC = new ArrayList<TrainingSequence<Character>>();
		List<TrainingSequence<MultipleAlignmentColumn>> LTSMA = new ArrayList<TrainingSequence<MultipleAlignmentColumn>>();		
		
		for (int j=0; j<data.size(); j++) {
			LTSC.add( data.get(j).getTrainingComponent("ref") );
			LTSMA.add( data.get(j).getTrainingComponent("aln") );			
		}
		
		
		predictorlp.train(LTSC);
		mo.train(model,LTSMA);
		
		for (int i=0; i<nFeatures; i++) {
			float[][] A = new float[span[i]][h.range()];
			logprob.add(A);
		}
		
		// Initialize an array to hold the Feature values which will be passed back:
		int nStates = model.getNumStates();
		transitions = new DenseBooleanMatrix2D[nFeatures];
		for (int i=0; i<nFeatures; i++) {
			transitions[i] = new DenseBooleanMatrix2D(nStates, nStates);
			for (int k=2; k<geometry.get(i).length; k+=2) {
				transitions[i].setQuick(geometry.get(i)[k], geometry.get(i)[k+1], true);
			}
		}
	
		for (int i=0; i<nFeatures; i++) {
			for (int j=0; j<span[i]; j++) {
				for (int k=0; k<h.range(); k++) {
					logprob.get(i)[j][k]=(float) 1.0;
				}
			}	
		}

		// In English, what I want to do is this.  Loop through all of the training data, once for each Feature.
		// While so doing, look for any positions where one of the allowed transitions for that feature occurs.
		// At such positions, increment the counts for logprob.
		for(TrainingSequence<Character> seq : LTSC) {
			int len = seq.length();
		
			for (int i=0; i<nFeatures; i++) {
				for (int ix=0; ix<len; ix++) {
					if ((ix>=offset[i]) && (ix-offset[i]+span[i] <= seq.length()) && (ix > 0) ) {
						int yprev = seq.getY(ix-1);
						int y = seq.getY(ix);
						for (int j=0; j<nTrans[i]; j++) {
							if ( (yprev == geometry.get(i)[2+2*j]) && (y == geometry.get(i)[2+2*j+1]) ) {
								for (int pos=0; pos<span[i]; pos++) {
									char c = seq.getX(ix - offset[i] + pos);
									logprob.get(i)[pos][h.hash(c)] += 1.0;
								}
							}
						}
					}
				}
			}	
		}
		
		// Above we got counts; we wanted logarithms; here is also where we normalize:		
		for (int i=0; i<nFeatures; i++) {			
			for (int j=0; j<span[i]; j++) {
				float norm = (float) 0.0;
				for (int k=0; k<h.range(); k++) {
					norm += (float) logprob.get(i)[j][k];
				}
				Assert.a(norm>0);
				for (int k=0; k<h.range(); k++) {
					logprob.get(i)[j][k] = (float) (Math.log(logprob.get(i)[j][k]) - Math.log(norm));
				}
			}	
		}

	}
}

