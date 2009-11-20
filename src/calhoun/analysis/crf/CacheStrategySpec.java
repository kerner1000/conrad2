/**
 * 
 */
package calhoun.analysis.crf;

import java.util.ArrayList;

import calhoun.util.Assert;

public class CacheStrategySpec {
	public enum CacheStrategy { COMPOSITE, CONSTANT, SPARSE, UNSPECIFIED, DENSE_NODE_BOUNDARY, DENSE, LENGTHFUNCTION /*, INTERVAL, LENGTH, NO_CACHE */}

	public CacheStrategy  strategy;
	public Object         details;
	
	public CacheStrategySpec(CacheStrategy strategy, Object details) {
		this.strategy = strategy;
		this.details = details;
	}
	
	public CacheStrategySpec(CacheStrategy strategy) {
		this(strategy,null);
	}
	
	/** Used in cases where the feature will return a value at every edge and/or node. */
	public static class DenseCachingDetails {
	
		public int nTables;
		public int nEvals;
		public int[] potential;
		public int[] tableNum;
		public short[] featureIndex;
		
		public void check() {
			Assert.a(nTables >= 0);
			Assert.a(potential.length == nEvals);
			Assert.a(tableNum.length == nEvals);
			Assert.a(featureIndex.length == nEvals);
			
			for (int j=0; j<nEvals; j++) {
				Assert.a(potential[j] >=0);
				Assert.a(tableNum[j]>=0);
				Assert.a(tableNum[j] < nTables);
				Assert.a(featureIndex[j] >=0);
			}
		}
	}
	
	// Used when a FeatureManagerNodeBoundaries is an explicit length node feature and the value it returns is obtained by subtracting two cumulative sums.
	public static class DenseBoundaryCachingDetails {
		
		public int nTables;
		public ArrayList<DenseBoundaryEntry> entries;
		
		public DenseBoundaryCachingDetails (int nTables) {
			this.nTables = nTables;
			entries = new ArrayList<DenseBoundaryEntry>();
		}
		
		public void add(int potential, int tableNum, int featureIndex, int rightPad, int leftPad) {
			entries.add(new DenseBoundaryEntry(potential,tableNum,featureIndex,rightPad,leftPad));
		}
		
		public void check() {
			int maxTable = 0;
			for (int j=0; j<entries.size(); j++) {
				DenseBoundaryEntry dbe = entries.get(j);
				dbe.check();
				int tnum = dbe.tableNum;
				if (tnum > maxTable) {
					maxTable = tnum;
				}			
			}
			Assert.a(maxTable == (nTables-1));
		}
	}
	
	public static class DenseBoundaryEntry {
		public int potential;
		public int tableNum;
		public int featureIndex;
		public int rightPad;
		public int leftPad;
		
		public DenseBoundaryEntry(int potential, int tableNum, int featureIndex, int rightPad, int leftPad) {
			this.potential = potential;
			this.tableNum = tableNum;
			this.featureIndex = featureIndex;
			this.rightPad = rightPad;
			this.leftPad = leftPad;
			check();
		}

		private void check() {
			Assert.a(potential >= 0);
			Assert.a(tableNum >= 0);
			Assert.a(featureIndex >= 0);
			Assert.a(rightPad >= 0); 
			Assert.a(leftPad >= 0);	
		}
	}
	
}