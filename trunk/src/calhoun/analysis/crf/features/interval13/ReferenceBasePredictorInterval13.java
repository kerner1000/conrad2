package calhoun.analysis.crf.features.interval13;

import calhoun.analysis.crf.CacheStrategySpec;
import calhoun.analysis.crf.FeatureManagerNodeBoundaries;
import calhoun.analysis.crf.CacheStrategySpec.CacheStrategy;


public class ReferenceBasePredictorInterval13 extends ReferenceBasePredictorInterval13Base implements FeatureManagerNodeBoundaries<Character> {
	
	private static final long serialVersionUID = -8460452348450096338L;

	public ReferenceBasePredictorInterval13() {
	}

	@Override
	public CacheStrategySpec getCacheStrategy() {
		
		CacheStrategySpec css = new CacheStrategySpec(CacheStrategy.DENSE_NODE_BOUNDARY);
		CacheStrategySpec.DenseBoundaryCachingDetails details = new CacheStrategySpec.DenseBoundaryCachingDetails(9); // we will use 9 tables
		

		// If you want to predict for real and mesh with a PWM, then set pads, for example, as below:
		//	donor[j]    = new PWMLookup(3,6,pseudoCounts);   // donor signal           xxx|GTxxxx 
		//	acceptor[j] = new PWMLookup(9,6,pseudoCounts);   // acceptor signal  xxxxxxxAG|xxxxxx
		//  start = new PWMLookup(9,6,pseudoCounts);             // start signal xxxxxxxxx|ATGxxx
		//  stop  = new PWMLookup(3,9,pseudoCounts);             // stop signal        xxx|TAGxxxxxx
		
		// IDEALLY HAVE SOME ASSERTIONS TYING THESE PADS TO THE PWM SIZES
		details.add(0,0,startIx + (multipleFeatures ? 0 : 0),Interval13Model.getPadIntergenic(),Interval13Model.getPadIntergenic());  // minlength 18
		details.add(1,1,startIx + (multipleFeatures ? 1 : 0),Interval13Model.getPadExon3prime(),Interval13Model.getPadExon5prime());  // minlength 9
		details.add(2,2,startIx + (multipleFeatures ? 1 : 0),Interval13Model.getPadExon3prime(),Interval13Model.getPadExon5prime());
		details.add(3,3,startIx + (multipleFeatures ? 1 : 0),Interval13Model.getPadExon3prime(),Interval13Model.getPadExon5prime());
		details.add(4,4,startIx + (multipleFeatures ? 2 : 0),Interval13Model.getPadIntron3prime(),Interval13Model.getPadIntron5prime());  // minlength 15
		details.add(5,4,startIx + (multipleFeatures ? 2 : 0),Interval13Model.getPadIntron3prime(),Interval13Model.getPadIntron5prime());
		details.add(6,4,startIx + (multipleFeatures ? 2 : 0),Interval13Model.getPadIntron3prime(),Interval13Model.getPadIntron5prime());
		details.add(7,5,startIx + (multipleFeatures ? 3 : 0),Interval13Model.getPadExon5prime(),Interval13Model.getPadExon3prime());  // minlength 9
		details.add(8,6,startIx + (multipleFeatures ? 3 : 0),Interval13Model.getPadExon5prime(),Interval13Model.getPadExon3prime());
		details.add(9,7,startIx + (multipleFeatures ? 3 : 0),Interval13Model.getPadExon5prime(),Interval13Model.getPadExon3prime());
		details.add(10,8,startIx + (multipleFeatures ? 4 : 0),Interval13Model.getPadIntron5prime(),Interval13Model.getPadIntron3prime());  // minlength 15
		details.add(11,8,startIx + (multipleFeatures ? 4 : 0),Interval13Model.getPadIntron5prime(),Interval13Model.getPadIntron3prime());
		details.add(12,8,startIx + (multipleFeatures ? 4 : 0),Interval13Model.getPadIntron5prime(),Interval13Model.getPadIntron3prime());		
		
		// NOTE: The minimum lengths of the Semi-Markov training and inference for each state must be AT least as big as the sum of the two pads at either end.
		
		details.check();
				
		css.details = details;
		
		return css;
	}
}
