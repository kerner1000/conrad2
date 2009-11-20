package flanagan.math.test;

import junit.framework.TestCase;
import flanagan.math.Minimisation;
import flanagan.math.MinimisationFunction;

public class NelderMeadTest extends TestCase {

	class Quadratic implements MinimisationFunction {
		int a = 2;
		int b = 3;
		int c = 4;
		public double function(double[] param) {
			// y = a(x-b)^2 + c
			// Minimum is (b,c)
			return a*Math.pow(param[0]-b,2)+c;
		}
	}
	
	class MultivariateQuadratic implements MinimisationFunction {
		int a = 2;
		int b = 3;
		int c = 4;
		int d = 5;
		public double function(double[] param) {
			// z = a(x-b)^2 + a(y-c)^2 + d
			// Minimum is (b,c,d)
			return a*Math.pow(param[0]-b,2)+a*Math.pow(param[1]-c,2)+d;
		}
	}
	
	public void testNelderMead() throws Exception {
		Minimisation m = new Minimisation();
		m.nelderMead(new Quadratic(), new double[] {1});
		assertEquals(4.0, m.getMinimum(), 0.0001);
		assertEquals(3.0, m.getParamValues()[0], 0.0001);

		m.nelderMead(new MultivariateQuadratic(), new double[] {1,1});
		assertEquals(5.0, m.getMinimum(), 0.0001);
		assertEquals(3.0, m.getParamValues()[0], 0.0001);
		assertEquals(4.0, m.getParamValues()[1], 0.0001);
	}

}
