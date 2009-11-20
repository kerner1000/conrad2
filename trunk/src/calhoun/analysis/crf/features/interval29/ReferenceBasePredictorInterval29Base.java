package calhoun.analysis.crf.features.interval29;

import java.util.List;

import calhoun.analysis.crf.AbstractFeatureManager;
import calhoun.analysis.crf.FeatureList;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.features.supporting.LogProbLookup;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.util.Assert;


public abstract class ReferenceBasePredictorInterval29Base extends AbstractFeatureManager<Character> {

	private static final long serialVersionUID = 8194502006226691957L;
	ModelManager model;
	int startIx;
		
	boolean multipleFeatures = false;
	
	double pseudoCounts;
	int lookback;
	
	LogProbLookup   intron;
	LogProbLookup   intergenic;
	LogProbLookup[] exonic;
	
	public ReferenceBasePredictorInterval29Base() {
	}

	public String getFeatureName(int featureIndex) {
		if(multipleFeatures) {
			int feat = featureIndex - startIx;
			String table = "";
			switch(feat) {
				case 0:
					table = "intergenic";
					break;
				case 1:
					table = "exon";
					break;
				case 2:
					table = "intron";
					break;
				case 3:
					table = "minus exon";
					break;
				case 4:
					table = "minus intron";
					break;
				default:
					Assert.a(false);
				}	
			return "referenceBasePredictorInterval29 "+table;
		}
		else {
			return "referenceBasePredictorInterval29";
		}
	}

	public int getNumFeatures() {
		return multipleFeatures ? 5 : 1;
	}
	
	public void train(int startingIndex, ModelManager modelInfo, List<? extends TrainingSequence<? extends Character>> data) {
		startIx = startingIndex;
		
		model = modelInfo;
		Interval29Tools.verify(modelInfo);

		pseudoCounts = 1.0;	
		lookback = 3;

		
		// Construct the space for the lookup tables.
		exonic = new LogProbLookup[3];
		for (int j=0; j<3; j++) {
			exonic[j] = new LogProbLookup(lookback,pseudoCounts);
		}
		intron     = new LogProbLookup(lookback,pseudoCounts);
		intergenic = new LogProbLookup(lookback,pseudoCounts);


		for(TrainingSequence<? extends Character> seq : data) {
			for (int pos=0; pos<seq.length(); pos++) {
				
				int state = seq.getY(pos);
				switch(state) {
				case(0):
				case(13):
				case(14):	
				case(21):
				case(22):					
					intergenic.increment(seq,pos,true);
					intergenic.increment(seq,pos,false);
					break;
				case(1):
				case(2):
				case(3):
					exonic[((pos-state+1)%3+3)%3].increment(seq,pos,true);
					break;
				case(4):
				case(5):
				case(6):
				case(15):
				case(16):
				case(17):
				case(18):
				case(19):
				case(20):					
					intron.increment(seq,pos,true);
					break;
				case(7):
				case(8):
				case(9):
					exonic[((-pos+state+1)%3+3)%3].increment(seq,pos,false);
					break;
				case(10):
				case(11):
				case(12):
				case(23):
				case(24):
				case(25):
				case(26):
				case(27):
				case(28):					
					intron.increment(seq,pos,false);
					break;
				default:
					Assert.a(false);
				}	
			}
		}
		
		for (int j=0; j<3; j++) {
			exonic[j].finalize();
		}
		intron.finalize();
		intergenic.finalize();
	}

	
	public void evaluateNode(InputSequence<? extends Character> seq, int pos, int state, FeatureList result) {
		double evaluation=0;

		int indexOffset = Integer.MIN_VALUE;
		int phase;
		switch(state) {
		case(0):
		case(13):
		case(14):	
		case(21):
		case(22):			
			evaluation = intergenic.lookup(seq,pos,true);
			indexOffset = 0;
			break;
		case(1):
		case(2):
		case(3):
			phase = ((pos-state+1)%3+3)%3;
			evaluation = exonic[phase].lookup(seq,pos,true);
			indexOffset = 1;// + phase;
			break;
		case(4):
		case(5):
		case(6):
		case(15):
		case(16):
		case(17):
		case(18):
		case(19):
		case(20):			
			evaluation = intron.lookup(seq,pos,true);
			indexOffset = 2;
			break;
		case(7):
		case(8):
		case(9):
			phase = ((-pos+state+1)%3+3)%3;
			evaluation = exonic[phase].lookup(seq,pos,false);
			indexOffset = 3;// + phase;
			break;
		case(10):
		case(11):
		case(12):
		case(23):
		case(24):
		case(25):
		case(26):
		case(27):
		case(28):			
			evaluation = intron.lookup(seq,pos,false);
			indexOffset = 4;
			break;		
		default:
			Assert.a(false);
		}
		
		result.addFeature(startIx + (multipleFeatures ? indexOffset : 0), evaluation);		
	}

	/** if true, a separate feature index is used for each state, creating 21 weights instead of 1.
	 * @return returns true if a separate feature index is used for each state
	 */
	public boolean isMultipleFeatures() {
		return multipleFeatures;
	}

	/**
	 * @param multipleFeatures The multipleFeatures to set.
	 */
	public void setMultipleFeatures(boolean weightPerState) {
		this.multipleFeatures = weightPerState;
	}

	
}
