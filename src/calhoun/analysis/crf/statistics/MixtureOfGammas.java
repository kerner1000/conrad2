package calhoun.analysis.crf.statistics;

import java.io.PrintStream;
import java.io.Serializable;
import java.util.Arrays;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.util.Assert;
import calhoun.util.ColtUtil;

public class MixtureOfGammas implements Serializable {
	private static final long serialVersionUID = 6269582005717276381L;
	private static final Log log = LogFactory.getLog(MixtureOfGammas.class);

	private double pdist1;
	private double shape1;
	private double lambda1;
	private double shape2;
	private double lambda2;

	public MixtureOfGammas(double pdist1, double shape1, double lambda1, double shape2, double lambda2) {
		Assert.a(pdist1>=0);
		Assert.a(pdist1<=1.0);
		Assert.a(shape1>0);
		Assert.a(lambda1>0);
		Assert.a(shape2>0);
		Assert.a(lambda2>0);
		
		this.pdist1 = pdist1;
		this.shape1 = shape1;
		this.lambda1 = lambda1;
		this.shape2 = shape2;
		this.lambda2 = lambda2;
	}


	
	public MixtureOfGammas(double[] lengths) {
		setup(lengths,false);
		log.info(summary());
	}

	public MixtureOfGammas(double[] lengths, boolean forceExponentialLength) {
		setup(lengths,forceExponentialLength);
		log.info(summary());
	}
	
	public void setup(double[] lengths, boolean forceExponentialLength) {
		int nLengths = lengths.length;		
		for (int j=0; j<nLengths; j++) {
			Assert.a(lengths[j]>0);
		}

		Assert.a(nLengths>=0);
		
		boolean exponentialDistribution = false;
		
		if (nLengths==0) {
			log.warn("Train mixture gamma called with no length inputs; returning an exponential distn with mean 100");
			double mean = 100;
			pdist1 = 1.0;
			shape1 = 1.0;
			lambda1 = 1/mean;
			shape2 = 1.0;
			lambda2 = 1.0;
			
			return;		
		}

		double[] sortedLengths = lengths;
		Arrays.sort(sortedLengths);
		if(log.isDebugEnabled())
			log.debug("Mixture of Gammas trainer called to model these lengths: "+ColtUtil.format(sortedLengths));
		
		
		if (forceExponentialLength) {
			log.warn("you called a mixture of gammas model but set flag to force it to model as an exponential length distribution.");
			exponentialDistribution = true;			
		}
		
		if (nLengths < 20) {
			log.warn("fewer than 20 lengths supplied for training; modeling with an exponential distribution instead of a mixture of Gammas");
			exponentialDistribution = true;
		}	
		
		boolean allSame = true;
		double firstLength = lengths[0];
		for (int j=0; j<lengths.length; j++) {
			if ((lengths[j]<0.99*firstLength) || (lengths[j]>1.01*firstLength)) {
				allSame = false;
			}
		}
		if (allSame) {
			log.warn("All the lengths we're asked to model are extremely close to the same value, suggesting the length is artifically fixed for some reason.  This is probably a mistake and some othe model of length is appropriate.  Will model as exponential length distribution.");
			exponentialDistribution = true;
		}
		
		if (exponentialDistribution) {
			double mean = BasicStats.meanDoubleArray(lengths);
			pdist1 = 1.0;
			shape1 = 1.0;
			lambda1 = 1/mean;
			shape2 = 1.0;
			lambda2 = 1.0;
			
			return;				
		}
		
		/* Think through in english and equations before trying to write out code
		 * 
		 * The hidden variables Y will be which of the two mixtures a particular point comes from.
		 * 
		 * The observed value X will be the observed intron length.
		 * 
		 * The parameters will be p (the probability of the biased coin saying distribution 1),
		 * the gamma parameters p1 and lambda1 for distribution 1, and the parameters p2 and
		 * lambda2 for other distribution.  Aggregating, let theta = (p,p1,lambda1,p2,lambda2).
		 * 
		 * The goal is to choose the parameters theta that maximize Pr(X|theta) = sum_Y ( Pr(Y,X|theta) ).
		 * 
		 * We do this using E-M algorithm, which iteratively increases Pr(X|theta) and will converge to
		 * a local but not necessarily global maximum.  Given theta_0 we need to be able to find
		 * theta_1 so that Pr(X|theta1) > Pr(X|theta0).
		 * 
		 * log[Pr(X|theta1)/Pr(X|theta0)] = log( sum_Y [ Pr(Y|X,theta0) * Pr(X,Y|theta1) / Pr(Y,X|theta0) ] )
		 *      [log is concave]    >=  sum_Y [ Pr(Y|X,theta0) * log [ Pr(X,Y|theta1) /  Pr(Y,X|theta0) ] ]
		 *                          = G(theta1) - G(theta0),
		 * where G(theta) = sum_Y [ Pr(Y|X,theta0) * log(Pr(X,Y|theta) ] is the auxillary function.
		 *      
		 * The EXPECTATION STEP is to form the function G in a way that can be evaluated as a function of theta.
		 *      
		 * The MAXIMIZATION STEP is to choose the value theta1 which maximizes G.  It must then be the
		 *      case that G(theta1)-G(theta0) >= 0, so that Pr(X|theta1) >= Pr(X|theta0). 
		 * 
		 * Need to introduce fake data, two points per distribution, to prevent ML convergence to prevent
		 * having one of the distributions converge on a single point, which is the non-useful trivial
		 * solution to the ML problem.  The fake data will be introduced as follows.
		 * We sort the observed lengths, and take the median.  Then the fake data for distribution 1
		 * is 0.90*median and 1.00*median; the fake data for distn 2 is 1.00*median and 1.10*median. 
		 */		

		// Introduce the fake data to prevent convergence to a spike; essentially a prior
		double median = BasicStats.medianDoubleArray(lengths);
		int len = nLengths+4;
		double[] x = new double[len];
		double[] post = new double[len];  // posterior probability of coming from distribution 1 not 2
		for (int j=0; j<nLengths; j++) { x[j] = lengths[j]; }
		x[nLengths]   = 0.90*median; post[nLengths]   = 1;  // fake data for dist 1
		x[nLengths+1] = 1.00*median; post[nLengths+1] = 1;  // fake data for dist 1
		x[nLengths+2] = 1.00*median; post[nLengths+2] = 0;  // fake data for dist 2
		x[nLengths+3] = 1.10*median; post[nLengths+3] = 0;  // fake data for dist 2
					
		
		// Start with initial seed guesses for p,p1,lambda1,p2,lambda2
		pdist1 = 0.5;
		shape1 = 15;
		lambda1 = 0.25;
		shape2 = 5;
		lambda2 = 0.05;
		
		
		for (int iteration=0; iteration<40; iteration++) {
		
			//System.out.println("Iteration " + iteration + "\tp=" + pdist1 + "\tshape1=" + shape1 + "\tlambda1=" + lambda1 + "\tmean1=" + shape1/lambda1 + "\tshape2=" + shape2 + "\tlambda2=" + lambda2 + "\tmean2=" + shape2/lambda2);
			
			
			// Do a single EM loop.
			// First the E-step:
			for (int j=0; j<nLengths; j++) {
				double p1 = pdist1*GammaDistribution.gamma(shape1,lambda1,x[j]);
				double p2 = (1-pdist1)*GammaDistribution.gamma(shape2,lambda2,x[j]);
				if (!(p1 + p2 > 0)) {
					Assert.a(false,"x[j]=" + x[j] + "  p1="+p1 + "  p2="+p2+"  pdist1=" + pdist1 + "  shape1="+shape1+"  lambda1="+ lambda1+"  shape2="+shape2+"  lambda2="+lambda2);
				}
				post[j] = p1 / (p1 + p2);
			}
			
			// And now the M-step.
			
			double mean1 = 0, mean2 = 0;
			double meanlog1 = 0, meanlog2 = 0;
			
			for (int j = 0; j < len; j++) {
				Assert.a(x[j] > 0);
				mean1 += post[j] * x[j];
				mean2 += (1 - post[j]) * x[j];
				meanlog1 += post[j] * Math.log(x[j]);
				meanlog2 += (1 - post[j]) * Math.log(x[j]);
			}
			pdist1 = BasicStats.meanDoubleArray(post);
			mean1 /= pdist1 * len;
			meanlog1 /= pdist1 * len;
			mean2 /= (1 - pdist1) * len;
			meanlog2 /= (1 - pdist1) * len;
			
			double[] plam1 = GammaDistribution.mleg(mean1, meanlog1);
			double[] plam2 = GammaDistribution.mleg(mean2, meanlog2);
			shape1 = plam1[0];
			lambda1 = plam1[1];
			shape2 = plam2[0];
			lambda2 = plam2[1];
			
		}
		
	}


	public double logEvaluate(double x) {
		double lret1 = GammaDistribution.lgamma(shape1,lambda1,x); // plus log(pdist1)
		double lret2 = GammaDistribution.lgamma(shape2,lambda2,x); // plus log(1-pdist1)
		
		if (pdist1>0.999) {
			return lret1;
		}
		if (pdist1<0.001) {
			return lret2;
		}
		
		double maxlog = Math.max(lret1,lret2);
		
		lret1 -= maxlog;
		lret2 -= maxlog;
		
		if (lret1<-100) {
			return maxlog + Math.log(1-pdist1) + lret2;
		}

		if (lret2<-100) {
			return maxlog + Math.log(pdist1) + lret1;
		}		
		
		double ret = maxlog + Math.log((pdist1)*Math.exp(lret1) + (1-pdist1)*Math.exp(lret2));	
		
		Assert.a((ret != Double.NEGATIVE_INFINITY) && (ret != Double.POSITIVE_INFINITY) && (!Double.isNaN(ret)));
		return ret;
	}
	

	public double evaluate(double x) {
		double ret = pdist1*(GammaDistribution.gamma(shape1,lambda1,x));
		ret += (1-pdist1)*(GammaDistribution.gamma(shape2,lambda2,x));
		Assert.a((ret != Double.NEGATIVE_INFINITY) && (ret != Double.POSITIVE_INFINITY) && (!Double.isNaN(ret)));
		return ret;
	}


	public void summarize(PrintStream out) {
		out.println(summary());
	
	}

	private String summary() {
		String ret = "";
		ret = ret + "MIXTURE OF GAMMAS INFO: pr(dist1)=" + pdist1;
		ret = ret + "  shape1=" + shape1;
		ret = ret + "  rate1=" + lambda1;
		ret = ret + "  mean1=" + (shape1/lambda1);
		ret = ret + "  shape2=" + shape2;
		ret = ret + "  rate2=" + lambda2;
		ret = ret + "  mean2=" + (shape2/lambda2);			
		return ret;
	}
	

	public double getMix() {
		return pdist1;
	}


	public double getMean() {
		double ret = pdist1*shape1/lambda1 + (1-pdist1)*shape2/lambda2;
		return ret;
	}
	
	
	
}


