package calhoun.analysis.crf.features.interval29;

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

public class PWMInterval29 extends AbstractFeatureManager<Character> implements FeatureManagerEdge<Character> {
	private static final long serialVersionUID = -7659288739348604129L;
	private static final Log log = LogFactory.getLog(PWMInterval29.class);
	boolean debug = log.isDebugEnabled();
	
	boolean multipleFeatures = false;
	
	int startIx;  // The index of the first feature managed by this FeatureManager
	ModelManager model;
	
	
	PWMLookup start,stop;
	PWMLookup[] donor,acceptor;
	double pseudoCounts;	

	
	public PWMInterval29() {
	}
	
	public int getNumFeatures() {
		return multipleFeatures ? 8 : 1;
	}	
	
	public String getFeatureName(int featureIndex) {
		if(multipleFeatures) {
			String feature = "";
			switch(featureIndex - startIx) {
			case 0:
				feature = "start";
				break;
			case 1:
				feature = "stop";
				break;
			case 2:
			case 3:
			case 4:
				feature = "donor"+(featureIndex - startIx - 2);
			case 5:
			case 6:
			case 7:
				feature = "acceptor"+(featureIndex - startIx - 5);
			}
			return "PwmFeatureInterval29 - "+feature;
		}
		else {
			Assert.a(featureIndex == startIx);
			return "PwmFeatureInterval29";
		}
	}

	public void train(int startingIndex, ModelManager modelInfo, List<? extends TrainingSequence<? extends Character>> data) {
		startIx = startingIndex;
		model = modelInfo;		
		
		Interval29Tools.verify(modelInfo);
				
		
		// Construct the space for the lookup tables.
		pseudoCounts = 1.0;	
		donor    = new PWMLookup[3];
		acceptor = new PWMLookup[3];
		for (int j=0; j<3; j++) {
			donor[j]    = new PWMLookup(Interval29Model.getPadExon3prime(),Interval29Model.getPadIntron5prime(),pseudoCounts);   // donor signal xxx|GTxxxx 
			acceptor[j] = new PWMLookup(Interval29Model.getPadIntron3prime(),Interval29Model.getPadExon5prime(),pseudoCounts);   // acceptor signal  xxxxxxxAG|xxxxxx
		}
		// Note: start PWM and stop PWM must extend equally far into the intergenic space, so that can set pads
		// stop and donor must also extend same amount into exon
		// start and acceptor must extend same amount into exon
		start = new PWMLookup(Interval29Model.getPadIntergenic(),Interval29Model.getPadExon5prime(),pseudoCounts);             // start signal xxxxxxxxx|ATGxxx
		stop  = new PWMLookup(Interval29Model.getPadExon3prime(),Interval29Model.getPadIntergenic(),pseudoCounts);             // stop signal xxx|TAGxxxxxx

				
		// Increment the lookup tables below
		for(TrainingSequence<? extends Character> seq : data) {
			for (int pos=1; pos<seq.length(); pos++) { // note start at one not zero, so can look back at prevState
				
				int state = seq.getY(pos);
				int prevState = seq.getY(pos-1);
				int iind;
				switch(Interval29Tools.edgeConstraints[prevState*Interval29Tools.numStates + state]) {
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
					// y(k) : k = i-j (mod 3)
					iind = Interval29Tools.check012(state-15);
					donor[iind].increment(seq,pos,true);
					break;
				case PACC:
					iind = Interval29Tools.check012(prevState-4);
					acceptor[iind].increment(seq,pos,true);
					break;
				case PSTOP:
					stop.increment(seq,pos,true);
					break;
				case MSTART:
					start.increment(seq,pos,false);
					break;
				case MDON:
					iind = Interval29Tools.check012(prevState-10);
					donor[iind].increment(seq,pos,false);
					break;
				case MACC:
					// y(k) : k = i+j (mod 3) 
					iind = Interval29Tools.check012(state-23);
					acceptor[iind].increment(seq,pos,false);
					break;
				case MSTOP:
					stop.increment(seq,pos,false);
					break;
				case PKEEPE:
				case PKEEPI:
				case MKEEPE:
				case MKEEPI:
				case PSTOPPED:
				case MSTARTED:
				case PWILLSTART:
				case MWILLSTOP:					
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
		switch(Interval29Tools.edgeConstraints[previousState*Interval29Tools.numStates + state]) {
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
			iind = Interval29Tools.check012(state-15);
			featureIndex = 2+iind;;
			val = donor[iind].lookup(seq,pos,true);
			break;
		case PACC:
			iind = Interval29Tools.check012(previousState-4);
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
			iind = Interval29Tools.check012(previousState-10);
			featureIndex = 2+iind;;
			val = donor[iind].lookup(seq,pos,false);
			break;
		case MACC:
			iind = Interval29Tools.check012(state-23);
			featureIndex = 5+iind;;
			val = acceptor[iind].lookup(seq,pos,false);
			break;
		case MSTOP:
			featureIndex = 1;
			val = stop.lookup(seq,pos,false);
			break;
		case PKEEPE:
		case PKEEPI:
		case MKEEPE:
		case MKEEPI:
		case PSTOPPED:
		case MSTARTED:
		case PWILLSTART:
		case MWILLSTOP:
			break;			
		default:
			Assert.a(false);  // We should have a complete enumeration of possibilities above.
		
		}
		Assert.a(val<=0);
		result.addFeature(startIx + (multipleFeatures ? featureIndex : 0),val);
	}

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

