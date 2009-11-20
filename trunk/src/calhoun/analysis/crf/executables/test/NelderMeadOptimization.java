package calhoun.analysis.crf.executables.test;

import java.util.Arrays;

import calhoun.util.ErrorException;
import flanagan.math.Minimisation;
import flanagan.math.MinimisationFunction;

public class NelderMeadOptimization {

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {	
				
		System.out.println("Example application of a Nelder-Meadfunction optimizer (without using gradient)");
		
		int maxIter = 5000;
		int nParm   = 10;
		
		MinimisationFunction mFunc = new MinimisationFunction() {
			public double function(double[] d) {
				double ret = 0;
				for (int j=0; j<d.length; j++) {
					ret += (d[j]-3)*(d[j]-3);
				}
				return ret;
			}
		};
		
		Minimisation m = new Minimisation();
		m.setNmax(maxIter);

		double[] starts = new double[nParm];
		Arrays.fill(starts, 1.0);

		double[] steps = new double[nParm];
		Arrays.fill(steps, 1.0);
	
		m.nelderMead(mFunc, starts, steps);
		if(!m.getConvStatus()) {
			throw new ErrorException("Convergence not reached.");
		}

		double[] results = m.getParamValues();
		
		System.out.println("Paremeters following optimization:\n");
		for (int j=0; j<nParm; j++) {
			System.out.println("  " + results[j]);
		}
		
		return;
	}

}
