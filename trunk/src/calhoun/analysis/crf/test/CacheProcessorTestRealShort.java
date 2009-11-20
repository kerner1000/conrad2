package calhoun.analysis.crf.test;

import calhoun.analysis.crf.Conrad;
import calhoun.util.AbstractTestCase;

// In this test want to use a 13tricycle model and train on
// real (but short) data, such as the ~10 genes in shortTrain.txt
// using a variety of cache processors, and make sure get identical
// answers (ie identical gradients or weights) in every case.

public class CacheProcessorTestRealShort extends AbstractTestCase {
	

	public void testEquivalentResultsVariousCacheProcessors() throws Exception {
		
		Conrad r = new Conrad("test/input/configFilesTricycle13/baselineTricycle13.xml");
		r.train("test/input/tinyTrain.tricycle13.txt");		

		Conrad rs = new Conrad("test/input/configFilesTricycle13/baselineTricycle13_smartcached.xml");
		rs.train("test/input/tinyTrain.tricycle13.txt");
		assertArrayEquals(rs.getWeights(), r.getWeights(), 0.0001);

		Conrad rsd = new Conrad("test/input/configFilesTricycle13/baselineTricycle13_smartcacheddense.xml");
		rsd.train("test/input/tinyTrain.tricycle13.txt");
		assertArrayEquals(rsd.getWeights(), r.getWeights(), 0.0001);
		
		Conrad rss = new Conrad("test/input/configFilesTricycle13/baselineTricycle13_smartcachedsparse.xml");
		rss.train("test/input/tinyTrain.tricycle13.txt");
		assertArrayEquals(rss.getWeights(), r.getWeights(), 0.0001);
		
		Conrad rc = new Conrad("test/input/configFilesTricycle13/baselineTricycle13_cached.xml");
		rc.train("test/input/tinyTrain.tricycle13.txt");
		assertArrayEquals(rc.getWeights(), r.getWeights(), 0.0001);
	}
	
	
	public void testCacheProcessorsWithInvalidation() throws Exception {
		
		Conrad rs = new Conrad("test/input/configFilesTricycle13/constrainingTricycle13_smartcached.xml");
		rs.train("test/input/tinyTrain.tricycle13.txt");
		
		Conrad rsd = new Conrad("test/input/configFilesTricycle13/constrainingTricycle13_smartcacheddense.xml");
		rsd.train("test/input/tinyTrain.tricycle13.txt");
		
		Conrad rss = new Conrad("test/input/configFilesTricycle13/constrainingTricycle13_smartcachedsparse.xml");
		rss.train("test/input/tinyTrain.tricycle13.txt");
		
		Conrad r = new Conrad("test/input/configFilesTricycle13/constrainingTricycle13.xml");
		r.train("test/input/tinyTrain.tricycle13.txt");
		
		Conrad rc = new Conrad("test/input/configFilesTricycle13/constrainingTricycle13_cached.xml");
		rc.train("test/input/tinyTrain.tricycle13.txt");
		
		assertArrayEquals(rc.getWeights(), r.getWeights(), 0.0001);
		assertArrayEquals(rs.getWeights(), r.getWeights(), 0.0001);
		assertArrayEquals(rsd.getWeights(), r.getWeights(), 0.0001);
		assertArrayEquals(rss.getWeights(), r.getWeights(), 0.0001);
	}
}
