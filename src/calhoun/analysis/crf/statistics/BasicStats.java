package calhoun.analysis.crf.statistics;

import java.util.Arrays;

import calhoun.util.Assert;

public class BasicStats {

	public static double meanDoubleArray( double[] x) {
		Assert.a(x.length>=1);
		double sum = 0;
		for (int j=0; j<x.length; j++) { sum+= x[j]; }
		return (sum/x.length);
	}

	public static double medianDoubleArray(double[] x) {
		Assert.a(x.length > 0);
		double[] y = x;
		Arrays.sort(y);
		return (y[y.length/2]);
	}

	public static double sumDoubleArray( double[] x) {
		double ret = 0.00;
		for (int j=0; j<x.length; j++) {
			ret += x[j];
		}
		return ret;
	}

	public static double L1Distance(double[] x, double[] y) {
		Assert.a(x.length == y.length);
		double ret = 0.00;
		for (int j=0; j<x.length; j++) {
			ret += Math.abs(x[j]-y[j]);
		}
		return ret;
	}

	public static int argmax(double[] y) {
		Assert.a(y.length>0);
		double val = y[0];
		int ret = 0;
		for (int j=1; j<y.length; j++) {
			if (y[j]>val)  {
				ret = j;
				val = y[j];
			}
		}
		return ret;
	}

	public static double max(double[] y) {
		Assert.a(y.length>0);
		double ret = y[0];
		for (int j=1; j<y.length; j++) {
			if (y[j]>ret)  {
				ret = y[j];
			}
		}
		return ret;
	}
	
}
