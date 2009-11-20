package calhoun.analysis.crf.test;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.statistics.BasicStats;
import calhoun.analysis.crf.statistics.GammaDistribution;
import calhoun.analysis.crf.statistics.MixtureOfGammas;
import calhoun.util.AbstractTestCase;
import calhoun.util.Assert;
import calhoun.util.FileUtil;

public class GammaDistributionsTest extends AbstractTestCase {
	private static final Log log = LogFactory.getLog(CRFIOTest.class);
	boolean debug = log.isDebugEnabled();

	
	
	public void testGammaDistribution() throws Exception {
		String fileName ="test/input/testIntronLengths.txt";
		double[] intronlengths = FileUtil.readDoublesFromSingleTabbedLine(fileName);
		int nIntrons = intronlengths.length;
		
		assertEquals(nIntrons,5863);
		assertEquals(intronlengths[0],59,0.00001);
		assertEquals(intronlengths[1],99,0.00001);
		
		double mean = 0.0;
		double meanlog = 0.0;
		for (int j=0; j<nIntrons; j++) {
			Assert.a(intronlengths[j] > 0);	
			mean += intronlengths[j];
			meanlog += Math.log(intronlengths[j]);
		}
		mean /= nIntrons;
		meanlog /= nIntrons;		
		System.out.println("Mean = " + mean);
		System.out.println("Meanlog = " + meanlog);
			
		double[] plam = GammaDistribution.mleg(mean,meanlog);
		double p = plam[0];
		double lam = plam[1];
		
		assertEquals(mean,p/lam);	
		System.out.println("Shape parameter p is "  +p);
		System.out.println("Rate parameter lambda is " + lam);
	}
	
	
	public void testGammaFunction() throws Exception {
		// for integers p, gamma(p) = (p-1)!
		assertEquals(GammaDistribution.gamma(1),1,0.000001);
		assertEquals(GammaDistribution.gamma(2),1,0.000001);
		assertEquals(GammaDistribution.gamma(3),2,0.000001);
		assertEquals(GammaDistribution.gamma(4),6,0.000001);
		assertEquals(GammaDistribution.gamma(5),24,0.000001);
		// Another identity is gamma(0.5) = sqrt(pi)
		assertEquals(GammaDistribution.gamma(0.5),Math.sqrt(Math.PI),0.000001);
	}
	

	
	
	public void testMixtureGammaModel() throws Exception {
		String fileName ="test/input/testIntronLengths.txt";
		double[] intronlengths = FileUtil.readDoublesFromSingleTabbedLine(fileName);
		int nIntrons = intronlengths.length;		
		
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
		 * solution to the ML problem.
		 */		

		// Introduce the fake data to prevent convergence to a spike; essentially a prior
		int len = nIntrons+4;
		double[] x = new double[len];
		double[] post = new double[len];  // posterior probability of coming from distribution 1 not 2
		for (int j=0; j<nIntrons; j++) { x[j] = intronlengths[j]; }
		x[nIntrons]   = 50; post[nIntrons]   = 1;  // fake data for dist 1
		x[nIntrons+1] = 60; post[nIntrons+1] = 1;  // fake data for dist 1
		x[nIntrons+2] = 60; post[nIntrons+2] = 0;  // fake data for dist 2
		x[nIntrons+3] = 90; post[nIntrons+3] = 0;  // fake data for dist 2
					
		
		// Start with initial guesses for p,p1,lambda1,p2,lambda2
		double p = 0.5;
		double shape1 = 15;
		double lambda1 = 0.25;
		double shape2 = 5;
		double lambda2 = 0.05;
		

		for (int iteration=0; iteration<20; iteration++) {
		
			System.out.println("Iteration " + iteration + "\tp=" + p + "\tshape1=" + shape1 + "\tlambda1=" + lambda1 + "\tmean1=" + shape1/lambda1 + "\tshape2=" + shape2 + "\tlambda2=" + lambda2 + "\tmean2=" + shape2/lambda2);
			
			
		// Do a single EM loop.
		// First the E-step:
		for (int j=0; j<nIntrons; j++) {
			double p1 = p*GammaDistribution.gamma(shape1,lambda1,x[j]);
			double p2 = (1-p)*GammaDistribution.gamma(shape2,lambda2,x[j]);
			Assert.a(p1+p2>0);
			post[j] = p1/(p1+p2); 
		}
		
		// And now the M-step.
		
		double mean1=0, mean2=0;
		double meanlog1=0, meanlog2=0;
		
		for (int j=0; j<len; j++) {
			Assert.a(x[j] > 0);	
			mean1 += post[j]*x[j];
			mean2 += (1-post[j])*x[j];
			meanlog1 += post[j]*Math.log(x[j]);
			meanlog2 += (1-post[j])*Math.log(x[j]);
		}
		p = BasicStats.meanDoubleArray(post);
		mean1 /= p*len;
		meanlog1 /= p*len;		
		mean2 /= (1-p)*len;
		meanlog2 /= (1-p)*len;
			
			
		double[] plam1 = GammaDistribution.mleg(mean1,meanlog1);
		double[] plam2 = GammaDistribution.mleg(mean2,meanlog2);
		shape1 = plam1[0];
		lambda1 = plam1[1];
		shape2 = plam2[0];
		lambda2 = plam2[1];
		
		}	
	}
	
	public void testMixtureGammaModelForceExponential() throws Exception {
		String fileName ="test/input/testIntronLengths.txt";
		double[] intronlengths = FileUtil.readDoublesFromSingleTabbedLine(fileName);
		new MixtureOfGammas(intronlengths); // tests that this runs to completion
		
		MixtureOfGammas mg = new MixtureOfGammas(intronlengths,true); // tests that runs to completion and is exponential
		
		for (int j=2; j<200; j++) {
			System.out.println("j="+j+"   " + Math.log(mg.evaluate((double) j)/mg.evaluate((double) (j-1))));
		}
	
		assertEquals(Math.log(mg.evaluate(2)/mg.evaluate(1)),Math.log(mg.evaluate(102)/mg.evaluate(101)),0.0001);
	}
	
	
	public void testMixtureOfGammasClass() throws Exception {
		String fileName ="test/input/testIntronLengths.txt";
		double[] intronlengths = FileUtil.readDoublesFromSingleTabbedLine(fileName);
	
		MixtureOfGammas MG = new MixtureOfGammas(intronlengths);
		MG.summarize(System.out);
		assertEquals(MG.getMix(),0.859,0.01);
		
		MixtureOfGammas MG2 = new MixtureOfGammas(0.86,71,1.27,4.1,0.041);		
		assertEquals(MG2.evaluate(55.0),0.053107,0.01);
		assertEquals(MG2.logEvaluate(55.0),-2.935446,0.01);
		
		// MG3 is exponential decay with mean 1
		MixtureOfGammas MG3 = new MixtureOfGammas(1.0, 1.0, 1.0, 1.0, 2.0);				
		assertEquals(MG3.evaluate(1.0),Math.exp(-1.0),0.01);
		assertEquals(MG3.logEvaluate(1.0),-1.0,0.01);
		
		// MG4 is exponential decay with mean 0.5
		MixtureOfGammas MG4 = new MixtureOfGammas(0, 1.0, 1.0, 1.0, 2.0);		
		assertEquals(MG4.evaluate(1.0),2.0*Math.exp(-2.0),0.01);
		assertEquals(MG4.logEvaluate(1.0),Math.log(2.0)-2.0,0.01);
		
		// MG5 is equal mix of above two distributions
		MixtureOfGammas MG5 = new MixtureOfGammas(0.5,1.0,1.0,1.0,2.0);		
		assertEquals(MG5.evaluate(1.0), 0.5*Math.exp(-1.0) + 0.5*2.0*Math.exp(-2.0),0.01);
		assertEquals(MG5.logEvaluate(1.0), Math.log(0.5*Math.exp(-1.0) + 0.5*2.0*Math.exp(-2.0)),0.01);
		
		// MG6 should be an exponential distribution with mean 3.5
		MixtureOfGammas MG6 = new MixtureOfGammas(new double[]{1,2,3,4,5,6});
		MG6.summarize(System.out);
		assertEquals(MG6.getMean(),3.5,0.001);
	}

}
