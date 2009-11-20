package calhoun.analysis.crf.features.tricycle13;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.AbstractFeatureManager;
import calhoun.analysis.crf.CacheStrategySpec;
import calhoun.analysis.crf.FeatureList;
import calhoun.analysis.crf.FeatureManagerEdge;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.BeanModel.Node;
import calhoun.analysis.crf.CacheStrategySpec.CacheStrategy;
import calhoun.analysis.crf.features.supporting.MarkovPredictorLogprob;
import calhoun.analysis.crf.features.tricycle13.EmissionMarkovFeature.MarkovHistory;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.seq.KmerHasher;
import calhoun.util.Assert;
import calhoun.util.DenseBooleanMatrix2D;

public class PositionWeightMatrixFeatures extends AbstractFeatureManager<Character> implements FeatureManagerEdge<Character> {
	private static final long serialVersionUID = -7659288739348604129L;
	private static final Log log = LogFactory.getLog(PositionWeightMatrixFeatures.class);
	boolean debug = log.isDebugEnabled();
	
	// Position weight matrices are used to model transitions between two extensive features,
	// for example a donor model for the transition between a positive stranded exon and a
	// positive stranded intron of the appropriate phase.
	
	// To specify a PWM, you specify its GEOMETRY, ie its span and where that span begins
	// relative to the transition itself (the offset), and the two hidden states before and
	// after the transition.
	
	// The features returned are that transition times the log probability of the observed
	// sequence within the window/span being modeled.  This must then be trained using
	// a TraningSequence<Character>.
	
	// Optionally, one may wish to subtract a double-counting correction if the observed
	// sequence in the window would have been predicted by something else.  For this you 
	// must also specify the predictor that would have been used by default for each base
	// (depending on the hidden state), and the sequence of hidden states to which this
	// predictor would have been applied over the span of the PWM.
	
	int startIx;  // The index of the first feature managed by this FeatureManager
	ModelManager model;
	
	// Following block are things that depend explicitly on and are calculated directly
	//  from geometry, included only for convenience. 
	int nFeatures;
	int[] span;
	int[] offset;
	int[] nTrans;
	KmerHasher   h; // for a single letter
	DenseBooleanMatrix2D[] transitions;
	
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

	
	
	// Following block is something that is only meaningful
	// if you're going to subtract the doublecounting correction
	// If you use this correction, you are currently required that each
	// Feature describes exactly one transition.
	boolean dcflag;
	List<int[]> dcc;  // DoubleCounting Correction
	MarkovPredictorLogprob predictorlp;

	
	transient InputSequence<? extends Character> lastSeq;
	int lastPos;
	float[] vals;
	
	boolean tieFlag = false;
	
	
	int UVCount = 0;
	
	public static class Geometry implements Serializable {
		private static final long serialVersionUID = 4896358213027322167L;

		int size;
		int transition;
		Node prev;
		Node current;

		List<Node> overlapCorrections;

		public Node getCurrent() {
			return current;
		}
		public void setCurrent(Node current) {
			this.current = current;
		}
		public List<Node> getOverlapCorrections() {
			return overlapCorrections;
		}
		public void setOverlapCorrections(List<Node> overlapCorrections) {
			this.overlapCorrections = overlapCorrections;
		}
		public Node getPrev() {
			return prev;
		}
		public void setPrev(Node prev) {
			this.prev = prev;
		}
		public int getSize() {
			return size;
		}
		public void setSize(int size) {
			this.size = size;
		}
		public int getTransition() {
			return transition;
		}
		public void setTransition(int transition) {
			this.transition = transition;
		}
	}

	List<Geometry> pwmGeometry;
	MarkovHistory markovHistory;
	
	///////////////////////////////////// Class variables above, methods below //////////
	public PositionWeightMatrixFeatures() { }
	
	public PositionWeightMatrixFeatures(List<int[]> geometry, List<int[]> dccorrection, List<int[]> markovhistory) {
		setThingsUp(geometry,dccorrection,markovhistory);
	}

	public PositionWeightMatrixFeatures(List<int[]> geometry2, List<int[]> dccorrection, List<int[]> markovhistory, List<int[]> flags) {
		tieFlag = true;
		setThingsUp(geometry2,dccorrection,markovhistory);
	}

	public void init() {
		List<int[]> geometry1 = new ArrayList(pwmGeometry.size());
		List<int[]> dccorrection = new ArrayList(pwmGeometry.size());
		for(Geometry g : pwmGeometry) {
			int[] params = new int[4];
			params[0] = g.getSize();
			params[1] = g.getTransition();
			params[2] = g.getPrev().getIndex();
			params[3] = g.getCurrent().getIndex();
			geometry1.add(params);
			int[] correction = new int[g.overlapCorrections.size()];
			for(int i=0; i<correction.length; ++i) {
				correction[i] = g.overlapCorrections.get(i).getIndex();
			}
			dccorrection.add(correction);
		}
		setThingsUp(geometry1, dccorrection, markovHistory.convert());
	}
	
	private void setThingsUp(List<int[]> geometry2, List<int[]> dccorrection, List<int[]> markovhistory) {
		this.predictorlp = new MarkovPredictorLogprob(markovhistory);
		
		setupGeometry(geometry2);	
		
		Assert.a(geometry.size()==nFeatures);
		Assert.a(dccorrection.size()==nFeatures);		
		for (int j=0; j<nFeatures; j++) {
			Assert.a(nTrans[j]==1);
			Assert.a(dccorrection.get(j).length == span[j]);
		}
		setupDoubleCountCorrections(dccorrection,predictorlp);		
	}
	
	
	private void setupDoubleCountCorrections(List<int[]> dccorrection, MarkovPredictorLogprob predictorlp) {
		dcflag = true;
		
		this.predictorlp = predictorlp;
		this.dcc = dccorrection;
	}

	private void setupGeometry(List<int[]> geometry) {
		this.geometry = geometry;
		nFeatures = geometry.size();
		span = new int[nFeatures];
		offset = new int[nFeatures];
		nTrans = new int[nFeatures];
		vals = new float[nFeatures];
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
	}

	public int getNumFeatures() {
		if (tieFlag) { return 1; }
		return nFeatures;
	}	
	
	public String getFeatureName(int featureIndex) {
		if (tieFlag) { return "tiedPwmFeature"; }
		
		int raw = featureIndex - startIx;
		int[] X = geometry.get(raw);
		
		String ret = "PWM.span" + X[0] + ".offset" + X[1];
		for (int j=2; j<X.length; j+=2) {
			ret = ret + ".(" + model.getStateName(X[j]) + "," + model.getStateName(X[j+1]) + ")"; 
		}	
		return ret;
	}


	
	//static int count = 0;
	public void evaluateEdge(InputSequence<? extends Character> seq, int pos, int previousState, int state, FeatureList result) {
		if(pos == 0) {
			return;
		}
		
		if((seq != lastSeq) || (pos != lastPos)) {
		//if((pos != lastPos)) {
			lastSeq = seq;
			lastPos = pos;
			updateVals(seq, pos);
			//System.out.println("Pos: "+pos+" Seq: "+seq);
			//if(++count > 10) {
			//	throw new RuntimeException();
			//}
		}
		
		if (tieFlag) {
			for (int j=0; j<nFeatures; j++) {
				if(transitions[j].getQuick(previousState, state)) {
					result.addFeature(startIx, vals[j]);
				}
			}		
		} else {
			for (int j=0; j<nFeatures; j++) {
				if(transitions[j].getQuick(previousState, state)) {
					result.addFeature(startIx + j, vals[j]);
				}
			}		
		}
	}

	void updateVals(InputSequence<? extends Character> seq, int ix) {
		UVCount++;
		for (int j=0; j<nFeatures; j++) {
			// Figure out what needed for Feature j, which might have several valid y pairs, at position i
			int[] geo = geometry.get(j);
			int spn = geo[0];
			int offset1 = geo[1];
			float val = 0;
			if ((ix>=offset1) && ((ix-offset1+spn)<=(seq.length())) ) {
				for (int i=0; i<spn; i++) {
					int pos = ix - offset1 + i;
					char c = seq.getX(pos);
					val = val + logprob.get(j)[i][h.hash(c)];
				}
			
				if (dcflag) {
					Assert.a(nTrans[j]==1);
					Assert.a(dcc.get(j).length == spn);
					
					for (int i=0; i<spn; i++) {
						val = val - predictorlp.logprob(dcc.get(j)[i],seq,ix-offset1+i);					
					}
				}
			}
			vals[j] = val;
		}			
	}

	public void train(int startingIndex, ModelManager modelInfo, List<? extends TrainingSequence<? extends Character>> data) {
		startIx = startingIndex;
		model = modelInfo;
	
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
		for(TrainingSequence<? extends Character> seq : data) {
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

		if (dcflag) {
			predictorlp.train(data);
		}
	}

	@Override
	public CacheStrategySpec getCacheStrategy() {
		return new CacheStrategySpec(CacheStrategy.SPARSE);
	}

	/**
	 * @return Returns the markovHistory.
	 */
	public MarkovHistory getMarkovHistory() {
		return markovHistory;
	}

	/**
	 * @param markovHistory The markovHistory to set.
	 */
	public void setMarkovHistory(MarkovHistory markovHistory) {
		this.markovHistory = markovHistory;
	}

	/**
	 * @return Returns the pwmGeometry.
	 */
	public List<Geometry> getPwmGeometry() {
		return pwmGeometry;
	}

	/**
	 * @param pwmGeometry The pwmGeometry to set.
	 */
	public void setPwmGeometry(List<Geometry> pwmGeometry) {
		this.pwmGeometry = pwmGeometry;
	}

	public boolean isTieFlag() {
		return tieFlag;
	}

	public void setTieFlag(boolean tieFlag) {
		this.tieFlag = tieFlag;
	}
}

