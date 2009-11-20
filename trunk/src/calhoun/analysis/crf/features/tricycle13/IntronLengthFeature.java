package calhoun.analysis.crf.features.tricycle13;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.AbstractFeatureManager;
import calhoun.analysis.crf.FeatureList;
import calhoun.analysis.crf.FeatureManagerNodeExplicitLength;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.analysis.crf.statistics.GammaDistribution;
import calhoun.util.Assert;

/**  Learns from data an explicit length distribution for introns, and also subtracts out
 * what it would have gotten from just the exponential distribution given by
 * DenseWeightedEdgeFeatures.
 */
public class IntronLengthFeature extends AbstractFeatureManager implements  FeatureManagerNodeExplicitLength {
	private static final long serialVersionUID = 8477631359065280630L;
	private static final Log log = LogFactory.getLog(IntronLengthFeature.class);
	boolean debug = log.isDebugEnabled();
	
	// Standard boring member variables
	int startIx;

	// Fundamental member variables
	double[] logProbExtend;       // given that in state i, what is log of probability of remaining in state i for another base?
	boolean[] explicitLengthFlag; // true if this is a state we wish to model explictly, otherwise false
	private boolean[] isIntron;	
		
	public String getFeatureName(int featureIndex) {
		Assert.a(featureIndex == startIx, "Invalid feature index");
		return "IntronLengthFeature";
	}

	public int getNumFeatures() { 
		return 1;
	}

	public void evaluateNodeLength(InputSequence seq, int pos, int length, int state, FeatureList result) {
		if(explicitLengthFlag[state]) {
			// Gaussian with mean mu=69 and stddev sigma=15 is (parms a guess from experience) is
			//  [1/(sigma*sqrt(2*pi))]  * exp[ - (x-mu)^2 / (2*sigma^2) ]		
			//double val = -3.6269887 - (length-69.0)*(length-69.0)/450.0 ;
			
			// Correct for what probably we would have gotten from the exponential state length.
			//val -= (length-1)*logProbExtend[state];				
			
			// Mixture of two gamma distributions with the following parameters:
			// p = 0.86 = probability from distribution 1:
			// dist1 is gamma distribution with shape 71 and lambda=1.27, so that mean=56
			// dist2 is gamma distribution with shape 4.1 and lambda=0.041, so that mean=101
			// [eventually this will be trained, but for now hardcoded.
			double pdist1 = 0.86;
			double shape1 = 71;
			double lambda1 = 1.27;
			double shape2 = 4.1;
			double lambda2 = 0.041;
			double val = pdist1*(GammaDistribution.gamma(shape1,lambda1,length));
			val += (1-pdist1)*(GammaDistribution.gamma(shape2,lambda2,length));
			val = Math.log(val);
			
			// Correct for what probably we would have gotten from the exponential state length.
			val -= (length-1)*logProbExtend[state];				

			result.addFeature(startIx, val);
		}
	}

	public void train(int startingIndex, ModelManager modelInfo, List data) {
		log.debug("Training the Feature for IntronLengths, explicitly modelled as a Gaussian");
		
		startIx = startingIndex;
		int nStates = modelInfo.getNumStates();
		
		// Say here which states modelled explicitly; the remainder are not
		explicitLengthFlag = new boolean[nStates];
		Arrays.fill(explicitLengthFlag, false);
		explicitLengthFlag[modelInfo.getStateIndex("intron1")]  = true;
		explicitLengthFlag[modelInfo.getStateIndex("intron2")]  = true;
		explicitLengthFlag[modelInfo.getStateIndex("intron3")]  = true;
		explicitLengthFlag[modelInfo.getStateIndex("intron1m")] = true;
		explicitLengthFlag[modelInfo.getStateIndex("intron2m")] = true;
		explicitLengthFlag[modelInfo.getStateIndex("intron3m")] = true;
		isIntron = explicitLengthFlag.clone();
		
		
		// Count transitions from the training data
		float[][] transitions = new float[nStates][nStates];
		for (int j=0; j<nStates; j++) {
			for (int k=0; k<nStates; k++) {
				transitions[j][k] = (float) 1.0; // pseudocounts
			}
		}
				
		//DoubleMatrix2D transitions  = new DenseDoubleMatrix2D(nStates, nStates);
		for(TrainingSequence seq : (List<TrainingSequence>) data) {
			// Start at 1 because there is no transition for the first element of the sequence.
			for(int pos = 1; pos < seq.length(); ++pos) { 
				int start = seq.getY(pos-1);
				int end = seq.getY(pos);
				transitions[start][end] += (float) 1.0; 
			}
		}

		logProbExtend = new double[nStates];
		for (int j=0; j<nStates; j++) {
			float rowtotal = (float) 0.0;
			for (int k=0; k<nStates; k++) {
				rowtotal += transitions[j][k];
			}
			logProbExtend[j] = (float) Math.log(transitions[j][j] / rowtotal);
		}
		
		log.debug("logprobextend for the variuos states are:");
		for (int j=0; j<modelInfo.getNumStates(); j++) {
			log.debug("  " + modelInfo.getStateName(j) + "   " + logProbExtend[j]);
		}
		
		ArrayList intronLengths = new ArrayList();
		//= new ArrayList<integer>;
		for(TrainingSequence seq : (List<TrainingSequence>) data) {
			int lastIntronStart = -1;
			int y = seq.getY(0);
			for(int pos = 1; pos < seq.length(); ++pos) { 
				int yprev = y;
				y = seq.getY(pos);
				
				if (isIntron[yprev] && (!isIntron[y]) && (lastIntronStart>=-1)) {
					intronLengths.add(pos-lastIntronStart);
					lastIntronStart = -1;
				}
				if ((!isIntron[yprev]) && (isIntron[y])) {
					lastIntronStart = pos;
				}
			}
		}
		
		log.debug("The intron lengths are:");
		for (int j=0; j<intronLengths.size(); j++) {
			log.debug(intronLengths.get(j) + ",");
		}
	}
}
