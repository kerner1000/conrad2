package calhoun.analysis.crf.executables.test;

import java.util.Arrays;

public class SpeedTest {
	static class Lookback {
		public Lookback() {
			indices = new int[10];
			Arrays.fill(indices, 2);
			values = new float[10];
			Arrays.fill(values, 10.4f);
		}
		int lookback;
		int[] indices;
		float[] values;
	}
	
	public static void main(String[] args) {
		Lookback[][] potentialArray = new Lookback[100][500];
		for(int pot=0; pot<potentialArray.length; ++pot) {
			Lookback[] lbArray = potentialArray[pot];
			for(int i= 0; i<lbArray.length; ++i) {
				lbArray[i] = new Lookback();
			}
		}
		
		long start = System.currentTimeMillis();
		double sum = 0.0;
		for(int iter = 0; iter<100; ++iter) {
			int len1 = potentialArray.length;
			for(int pot=0; pot<len1; ++pot) {
				Lookback[] lbArray = potentialArray[pot];
				int len2 = lbArray.length;
				for(int i=0; i<len2; ++i) {
					int len3 = lbArray[i].indices.length;
					int[] lbIndices = lbArray[i].indices;
					float[] lbValues = lbArray[i].values;
					for(int j = 0; j<len3; ++j) {
						sum += lbIndices[j]*lbValues[j];
					}
				}
			}
		}
		System.out.println(System.currentTimeMillis() - start);

		int[][] lengths = new int[100][500];
		int[][][] arIndices = new int[100][500][10];
		float[][][] arValues= new float[100][500][10];

		start = System.currentTimeMillis();
		for(int iter = 0; iter<100; ++iter) {
			sum = 0.0;
			int len1 = lengths.length;
			for(int pot=0; pot<len1; ++pot) {
				int[] lbArray = lengths[pot];
				int len2 = lbArray.length;
				for(int i=0; i<len2; ++i) {
					int len3 = arIndices[pot][i].length;
					int[] lbIndices = arIndices[pot][i];
					float[] lbValues = arValues[pot][i];
					for(int j = 0; j<len3; ++j) {
						sum += lbIndices[j]*lbValues[j];
					}
				}
			}
		}
		System.out.println(System.currentTimeMillis() - start);
	}

}
