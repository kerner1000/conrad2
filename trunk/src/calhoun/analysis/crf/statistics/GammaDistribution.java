package calhoun.analysis.crf.statistics;

import calhoun.util.Assert;

public class GammaDistribution {

	/* See http://mathworld.wolfram.com/GammaDistribution.html
	 * 
	 * 
	 * 
	 * 
	 */
	
	public static double lgamma(double p, double lambda, double x) {
		Assert.a(lambda>0);
		Assert.a(x>0);
		Assert.a(p>0);
		double ret = (p-1)*Math.log(x) -x*lambda + p*Math.log(lambda) - lgamma(p);
		if (Double.isNaN(ret)) {
			Assert.a(false,"  pow1=" +Math.pow(x,p-1)+"  exp1="+Math.exp(-x*lambda)+"  pow2="+Math.pow(lambda,p) + "  gamma="+gamma(p) );
		}
		return ret;
	}
	
	
	public static double gamma(double p, double lambda, double x) {
		Assert.a(lambda>0);
		Assert.a(x>0);
		Assert.a(p>0);

		double lret = lgamma(p,lambda,x);
		double ret = Math.exp(lret);
		//Assert.a(!Double.isNaN(ret));
		Assert.a((ret != Double.NEGATIVE_INFINITY) && (ret != Double.POSITIVE_INFINITY) && (!Double.isNaN(ret)));
		return ret;
	}

	public static double gamma(double x) {
		/* The gamma function, which is the improper integral
		 * integral_0^infinity t^(x-1) * exp(-t) dt
		 * 
		 * 	x*gamma(x) = gamma(x+1).
		 * 
		 * For integers x>=1, gamma(x) = (x-1)! 
		 */
		Assert.a(x>0);
		double y = lgamma(x);
		return Math.exp(y);
	}
	
	
	 private static double lgamma(double x) {
		 Assert.a(x>0);
		 
		 double l2pi = Math.log(2*Math.PI);
		 double term, y, zz ; 
		 int k, t ;

		 if (x<10.0) return (lgamma(x+1.0) - Math.log(x)) ;
		 y = (x-0.5)*Math.log(x) -x + 0.5*l2pi ;  
		 zz = 1.0/x ;
		 for (k=1; (2*k)<= bern.length ; k++)  {
		  t = 2*k ;
		  term = bern[2*k]/(double) (t*(t-1)) ;
		  term *= zz ;
		  y += term ;
		  zz /= (x*x) ;
		 }
		 
		return y;
	}


	 

	// Ported from Nick Patterson's code 20060203
	 static private double[] bern = new double[]{1.0,
			 -1.0/2.0,
			 1.0/6.0,
			 0,
			 -1.0/30.0,
			 0,
			 1.0/42.0,
			 0,
			 -1.0/30.0,
			 0,
			 5.0/66.0,
			 0,
			 -691.0/2730.0,
			 0,
			 7.0/6.0
	 };
		
	 // Ported from Nick Patterson's code 20060203
	 public static double[] mleg(double a1, double a2) // both p and lam assumed of size 1 
	//  solve 
	//  p/lam = a1 ; psi(p) - log(lam) = a2 ;
	//  Thus psi(p) - log(p) = a2 - log(a1) 
	{
	   int iter ;
	   double s, pp, ll  ;
	   double top, bot, fval ;

	
	  s = a2 - Math.log(a1) ;   

	  Assert.a(s<=0.0 , "log E(x) < E(log (x)) \n" ) ;
	  
	  pp = -s ;

	  for (iter = 1; iter <= 30 ; ++iter) {  
	   fval = s - (psi(pp) - Math.log (pp)) ;
	   if (fval<0.0)  break ;
	   pp *= 2.0 ;
	  }

	  for (iter = 1; iter <= 30 ; ++iter) {  
	   fval = s - (psi(pp) - Math.log (pp)) ;
	   if (fval>0.0)  break ;
	   pp /= 2.0 ;
	  }
	  
	  for (iter = 1; iter <= 10 ; ++iter) {  
	   fval = psi(pp) - Math.log (pp) ;
	   top = s-fval ;
	   bot =  tau(pp) - (1.0/pp) ;
	   pp += top/bot ;
	  }
	  ll = pp/a1 ;
	  
	  return (new double[]{pp,ll});
	}



	private static double psi(double x) {
		// Ported from Nick Pattersons code 20060203
		// psi is the derivative of the gamma function
		
		 double y, zz, term ;
		 int k ;

		 Assert.a( x > 0.0 , "(psi) bad value: " +  x) ;

		 Assert.a(bern.length == 15);
		
		 if (x<10.0) return (psi(x+1.0) - 1.0/x) ;

		 y = Math.log(x) - 1.0/(2.0*x) ;
		 zz = 1.0 ;
		 for (k=1; (2*k)< bern.length ; k++)  {
		  zz /= (x*x) ;
		  term = bern[2*k]/(double) (2*k) ;
		  term *= zz ;
		  y -= term ;
		 }
		 return y ;
	}


	private static double tau(double x)
	//	 Ported from Nick Pattersons code 20060203
	//	derivative of psi 	
	{
		 double y, zz, term ;
		 int k ;

		 Assert.a(x>0.0 , "(tau) bad value: " +  x) ;
		 
		 if (x<10.0) return (tau(x+1.0) + 1.0/(x*x)) ;

		 y = 1.0/x  + 1.0/(2.0*x*x) ;
		 zz = 1.0/x ;
		 for (k=1; (2*k)< bern.length ; k++)  {
		  zz /= (x*x) ;
		  term = bern[2*k]/(double) (2*k) ;
		  term *= zz ;
		  term *= - (double) (2*k) ;
		  y -= term ;
		 }
		 return y ;
	}


}
