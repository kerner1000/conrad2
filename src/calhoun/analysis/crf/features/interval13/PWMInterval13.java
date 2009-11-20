package calhoun.analysis.crf.features.interval13;

import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.AbstractFeatureManager;
import calhoun.analysis.crf.CacheStrategySpec;
import calhoun.analysis.crf.FeatureList;
import calhoun.analysis.crf.FeatureManagerEdge;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.CacheStrategySpec.CacheStrategy;
import calhoun.analysis.crf.features.supporting.PWMLookup;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.util.Assert;

public class PWMInterval13 extends AbstractFeatureManager<Character> implements FeatureManagerEdge<Character> {
	private static final long serialVersionUID = -7659288739348604129L;
	private static final Log log = LogFactory.getLog(PWMInterval13.class);
	boolean debug = log.isDebugEnabled();
	
	boolean multipleFeatures = false;
	
	int startIx;  // The index of the first feature managed by this FeatureManager
	ModelManager model;
	
	
	PWMLookup start,stop;
	PWMLookup[] donor,acceptor;
	double pseudoCounts;	

	
	public PWMInterval13() {
	}
	
	public int getNumFeatures() {
		return multipleFeatures ? 8 : 1;
	}	
	
	public String getFeatureName(int featureIndex) {
		if(multipleFeatures) {
			String[] vals = new String[] { "Start", "Stop", "Donor 0", "Donor 1", "Donor 2", "Acceptor 0", "Acceptor 1", "Acceptor 2"};
			String feature = vals[featureIndex - startIx];
			return feature + " PWM";
		}
		else {
			Assert.a(featureIndex == startIx);
			return "PwmFeatureInterval13";
		}
	}

	public void train(int startingIndex, ModelManager modelInfo, List<? extends TrainingSequence<? extends Character>> data) {
		startIx = startingIndex;
		model = modelInfo;		
		
		Interval13Tools.verify(modelInfo);
				
		
		// Construct the space for the lookup tables.
		pseudoCounts = 1.0;	
		donor    = new PWMLookup[3];
		acceptor = new PWMLookup[3];
		for (int j=0; j<3; j++) {
			donor[j]    = new PWMLookup(Interval13Model.getPadExon3prime(),Interval13Model.getPadIntron5prime(),pseudoCounts);   // donor signal xxx|GTxxxx 
			acceptor[j] = new PWMLookup(Interval13Model.getPadIntron3prime(),Interval13Model.getPadExon5prime(),pseudoCounts);   // acceptor signal  xxxxxxxAG|xxxxxx
		}
		// Note: start PWM and stop PWM must extend equally far into the intergenic space, so that can set pads
		// stop and donor must also extend same amount into exon
		// start and acceptor must extend same amount into exon
		start = new PWMLookup(Interval13Model.getPadIntergenic(),Interval13Model.getPadExon5prime(),pseudoCounts);             // start signal xxxxxxxxx|ATGxxx
		stop  = new PWMLookup(Interval13Model.getPadExon3prime(),Interval13Model.getPadIntergenic(),pseudoCounts);             // stop signal xxx|TAGxxxxxx

				
		// Increment the lookup tables below
		for(TrainingSequence<? extends Character> seq : data) {
			for (int pos=1; pos<seq.length(); pos++) { // note start at one not zero, so can look back at prevState
				
				int state = seq.getY(pos);
				int prevState = seq.getY(pos-1);
				int iind;
				switch(Interval13Tools.edgeConstraints[prevState*Interval13Tools.numStates + state]) {
				case NONE:
				case PCODE:
				case MCODE:
					break;
				case NEVER:
					Assert.a(false,"pos = "+pos+" prevState = " + modelInfo.getStateName(prevState) + "   State = " + modelInfo.getStateName(state));  // A nice side effect of making sure the input sequence is legal, can omit this if you want to.
					break;
				case PSTART:
					start.increment(seq, pos,true);
					break;
				case PDON:
					iind = Interval13Tools.check012(state-4);
					donor[iind].increment(seq,pos,true);
					break;
				case PACC:
					iind = Interval13Tools.check012(prevState-4);
					acceptor[iind].increment(seq,pos,true);
					break;
				case PSTOP:
					stop.increment(seq,pos,true);
					break;
				case MSTART:
					start.increment(seq,pos,false);
					break;
				case MDON:
					iind = Interval13Tools.check012(prevState-10);
					donor[iind].increment(seq,pos,false);
					break;
				case MACC:
					iind = Interval13Tools.check012(state-10);
					acceptor[iind].increment(seq,pos,false);
					break;
				case MSTOP:
					stop.increment(seq,pos,false);
					break;
				default:
					Assert.a(false);  // We should have a complete enumeration of possibilities above.
				}
			}
		}
		
		for (int j=0; j<3; j++) {
			donor[j].completeCounts();
			acceptor[j].completeCounts();
		}
		start.completeCounts();
		stop.completeCounts();
	}
	
	
	public void evaluateEdge(InputSequence<? extends Character> seq, int pos, int previousState, int state, FeatureList result) {
		double val = 0.0;	
		
		int featureIndex = Integer.MIN_VALUE;
		int iind;
		switch(Interval13Tools.edgeConstraints[previousState*Interval13Tools.numStates + state]) {
		case NONE:
		case PCODE:
		case MCODE:
			break;
		case NEVER:
			Assert.a(false);  // A nice side effect of making sure the input sequence is legal, can omit this if you want to.
			break;
		case PSTART:
			featureIndex = 0;
			val = start.lookup(seq, pos,true);
			break;
		case PDON:
			iind = Interval13Tools.check012(state-4);
			featureIndex = 2+iind;;
			val = donor[iind].lookup(seq,pos,true);
			break;
		case PACC:
			iind = Interval13Tools.check012(previousState-4);
			featureIndex = 5+iind;;
			val = acceptor[iind].lookup(seq,pos,true);
			break;
		case PSTOP:
			featureIndex = 1;
			val = stop.lookup(seq,pos,true);
			break;
		case MSTART:
			featureIndex = 0;
			val = start.lookup(seq,pos,false);
			break;
		case MDON:
			iind = Interval13Tools.check012(previousState-10);
			featureIndex = 2+iind;;
			val = donor[iind].lookup(seq,pos,false);
			break;
		case MACC:
			iind = Interval13Tools.check012(state-10);
			featureIndex = 5+iind;;
			val = acceptor[iind].lookup(seq,pos,false);
			break;
		case MSTOP:
			featureIndex = 1;
			val = stop.lookup(seq,pos,false);
			break;
		default:
			Assert.a(false);  // We should have a complete enumeration of possibilities above.
		
		}
		Assert.a(val<=0);
		result.addFeature(startIx + (multipleFeatures ? featureIndex : 0),val);
	}

	@Override
	public CacheStrategySpec getCacheStrategy() {
		return new CacheStrategySpec(CacheStrategy.UNSPECIFIED);
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

