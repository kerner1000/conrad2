package calhoun.analysis.crf.features.interval29;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.CacheStrategySpec;
import calhoun.analysis.crf.FeatureList;
import calhoun.analysis.crf.FeatureManagerNodeExplicitLength;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.CacheStrategySpec.CacheStrategy;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.analysis.crf.io.SequenceConverter;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.analysis.crf.statistics.MixtureOfGammas;
import calhoun.util.Assert;

public class StateLengthLogprobInterval29 implements FeatureManagerNodeExplicitLength<Character> {	
	private static final long serialVersionUID = 8685543199212865835L;

	/* Returns the log probability of a state having the specified duration
	 * returned value depends only on state and length
	 * 
	 * trained using mixture of gammas model, with a few tweaks:
	 *   a) if all values identical (eg always have exactly 200 intergenic bases in training examples), then model as exponential length
	 *   b) Otherwise introduce priors at 90 and 110% of median and at 95% and 105% of median, and train mixture of gammas using EM.
	 *      This prevents convergence to a case where one of the two components is a single data point and likelihood is infinite.
	 *     
	 * Some notes on normalization: (see also mixGamma):
	 *    If you have a given state and evaluate these probabilities (before taking logs) and add up over all possible lengths,
	 *    the result should be approximately one, but not exactly.  Two reasons for discrepancy:
	 *    a) mixGamma, integrated from to to infinity, should add up to 1.0.  But summing at 1,2,3,4,... is a discrete
	 *      approximation and might not agree exactly.
	 *    b) The summation isn't from 1,2,3,4,... but from minLength[state] to maxLength[state].  However, this feature does not
	 *      have access to that min/max length information, so normalization for this reason must happen downstream.
	 */
	private static final Log log = LogFactory.getLog(StateLengthLogprobInterval29.class);
	
	private int startIx;
	ModelManager mi;
	private String inputComponentName;

	MixtureOfGammas intergenicMixGamma;
	MixtureOfGammas exonMixGamma;
	MixtureOfGammas intronMixGamma;
	
	private boolean forceExponential = false;
	private boolean exonExponential = false;
	private boolean multipleFeatures = true;

	public void setForceExponential(boolean forceExponential) {
		this.forceExponential = forceExponential;
	}
	
	public void setExonExponential(boolean exonExponential) {
		this.exonExponential = exonExponential;
	}
	
	public void evaluateNodeLength(InputSequence<? extends Character> seq, int pos, int length, int state, FeatureList result) {
		Assert.a(length>0);
		MixtureOfGammas mg = null;
		int indexOffset = Integer.MIN_VALUE;
		switch (state) {
		case(0):
			indexOffset = 0;
			mg = intergenicMixGamma;
		break;
		case(1):
		case(2):
		case(3):
		case(7):
		case(8):
		case(9):
			indexOffset = 1;
			mg = exonMixGamma;
		break;
		case(4):
		case(5):
		case(6):
		case(10):
		case(11):
		case(12):
			indexOffset = 2;
			mg = intronMixGamma;
		break;
		case(13):
		case(14):
		case(15):
		case(16):
		case(17):
		case(18):
		case(19):
		case(20):
		case(21):
		case(22):
		case(23):
		case(24):
		case(25):
		case(26):
		case(27):
		case(28):			
			indexOffset = 3;
			break;
		default:
			Assert.a(false);
		}
		double val;
		if (mg != null) {
			val = mg.logEvaluate((double) length);
		} else {
			Assert.a(state >= 13 && state <= 28, "Trying to use exact state lengths for non-boundary state.");
			// XXX: Should we be doing this?
			if (length == 2) {
				val = -0.5;
			} else {
				val = -200;
			}
		}
		Assert.a((val != Double.NEGATIVE_INFINITY) && (val != Double.POSITIVE_INFINITY) && (!Double.isNaN(val)));
		if (val>0) {
			log.warn("About to return a state length logprob evaluation that is greater than zero, see notes in source code.");
		}
		result.addFeature(startIx + (multipleFeatures ? indexOffset : 0),val);
	}

	public CacheStrategySpec getCacheStrategy() {
		return new CacheStrategySpec(CacheStrategy.LENGTHFUNCTION);
	}

	public String getFeatureName(int featureIndex) {
		if(multipleFeatures) {
			String type = "";
			switch(featureIndex - startIx) {
			case 0:
				type = "intergenic";
				break;
			case 1:
				type = "exon";
				break;
			case 2:
				type = "intron";
				break;
			case 3:
				type = "boundary";
				break;				
			}
			return "StateDurationLogProbForModelInterval21 "+type;
		}
		else {
			return "StateDurationLogProbForModelInterval21";
		}
	}

	public String getInputComponent() {
		return inputComponentName;
	}

	public void setInputComponent(String name) {
		inputComponentName = name;		
	}

	public void train(int startingIndex, ModelManager modelInfo, List<? extends TrainingSequence<? extends Character>> data) {
		startIx = startingIndex;
		mi = modelInfo;
		Assert.a(mi.getNumStates()==29);
		
		ArrayList<ArrayList<Integer>> stateDurations;
		
		stateDurations = SequenceConverter.stateVector2StateLengths(data,mi.getNumStates());
		//stateDurations = SequenceConverter.addInInterval29Boundaries(stateDurations);
		List<Integer> exonLengths = new ArrayList<Integer>();
		List<Integer> intronLengths = new ArrayList<Integer>();
		List<Integer> intergenicLengths = new ArrayList<Integer>();
		
		intergenicLengths.addAll(stateDurations.get(0));
		
		exonLengths.addAll(stateDurations.get(1));
		exonLengths.addAll(stateDurations.get(2));
		exonLengths.addAll(stateDurations.get(3));
		exonLengths.addAll(stateDurations.get(7));
		exonLengths.addAll(stateDurations.get(8));
		exonLengths.addAll(stateDurations.get(9));
		
		intronLengths.addAll(stateDurations.get(4));
		intronLengths.addAll(stateDurations.get(5));
		intronLengths.addAll(stateDurations.get(6));
		intronLengths.addAll(stateDurations.get(10));
		intronLengths.addAll(stateDurations.get(11));
		intronLengths.addAll(stateDurations.get(12));
		
		double[] inter = new double[intergenicLengths.size()];
		for (int j=0; j<intergenicLengths.size(); j++) {
			inter[j] = (double) intergenicLengths.get(j);
		}
		
		double[] exon = new double[exonLengths.size()];
		for (int j=0; j<exonLengths.size(); j++) {
			exon[j] = (double) exonLengths.get(j);
		}
		
		double[] intron = new double[intronLengths.size()];
		for (int j=0; j<intronLengths.size(); j++) {
			intron[j] = (double) intronLengths.get(j);
		}

		if (forceExponential) {
			intergenicMixGamma = new MixtureOfGammas(inter,true);
			exonMixGamma       = new MixtureOfGammas(exon,true);
			intronMixGamma     = new MixtureOfGammas(intron,true);
		} else if (exonExponential) {
			intergenicMixGamma = new MixtureOfGammas(inter,true);  
			exonMixGamma       = new MixtureOfGammas(exon,true);
			intronMixGamma     = new MixtureOfGammas(intron);				
		} else {
			// by default, only intergenic regions modeled with exp length distributions
			intergenicMixGamma = new MixtureOfGammas(inter,true);  
			exonMixGamma       = new MixtureOfGammas(exon);
			intronMixGamma     = new MixtureOfGammas(intron);			
		}
	}

	public int getNumFeatures() {
		return multipleFeatures ? 4 : 1;
	}

	/** 
	 * @return Returns the multipleFeatures.
	 */
	public boolean isMultipleFeatures() {
		return multipleFeatures;
	}

	/** set to true to indicate that intergenic, exonic, and intergenic lengths should each get a separate weight.
	 * @param multipleFeatures The multipleFeatures to set.
	 */
	public void setMultipleFeatures(boolean multipleFeatures) {
		this.multipleFeatures = multipleFeatures;
	}
}
