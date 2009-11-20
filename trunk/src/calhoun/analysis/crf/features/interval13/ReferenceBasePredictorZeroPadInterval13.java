package calhoun.analysis.crf.features.interval13;

import calhoun.analysis.crf.CacheStrategySpec;
import calhoun.analysis.crf.FeatureManagerNodeBoundaries;
import calhoun.analysis.crf.CacheStrategySpec.CacheStrategy;


public class ReferenceBasePredictorZeroPadInterval13 extends ReferenceBasePredictorInterval13Base implements FeatureManagerNodeBoundaries<Character> {
	
	private static final long serialVersionUID = -1104446043334637758L;

	public ReferenceBasePredictorZeroPadInterval13() {
	}

	@Override
	public CacheStrategySpec getCacheStrategy() {
		
		CacheStrategySpec css = new CacheStrategySpec(CacheStrategy.DENSE_NODE_BOUNDARY);
		CacheStrategySpec.DenseBoundaryCachingDetails details = new CacheStrategySpec.DenseBoundaryCachingDetails(9); // we will use 9 tables
		
		// state, table, featureIx, rightPad, leftPad
		// If you want to debug and verify equality with the Markov case, then set all the pads to zero as below:
		
		details.add(0,0,startIx + (multipleFeatures ? 0 : 0),0,0);
		details.add(1,1,startIx + (multipleFeatures ? 1 : 0),0,0);
		details.add(2,2,startIx + (multipleFeatures ? 1 : 0),0,0);
		details.add(3,3,startIx + (multipleFeatures ? 1 : 0),0,0);
		details.add(4,4,startIx + (multipleFeatures ? 2 : 0),0,0);
		details.add(5,4,startIx + (multipleFeatures ? 2 : 0),0,0);
		details.add(6,4,startIx + (multipleFeatures ? 2 : 0),0,0);
		details.add(7,5,startIx + (multipleFeatures ? 3 : 0),0,0);
		details.add(8,6,startIx + (multipleFeatures ? 3 : 0),0,0);
		details.add(9,7,startIx + (multipleFeatures ? 3 : 0),0,0);
		details.add(10,8,startIx + (multipleFeatures ? 4 : 0),0,0);
		details.add(11,8,startIx + (multipleFeatures ? 4 : 0),0,0);
		details.add(12,8,startIx + (multipleFeatures ? 4 : 0),0,0);

		details.check();
				
		css.details = details;
		
		return css;
	}
}
