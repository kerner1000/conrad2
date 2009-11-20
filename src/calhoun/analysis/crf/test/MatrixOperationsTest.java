package calhoun.analysis.crf.test;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.util.AbstractTestCase;
import calhoun.util.ColtUtil;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

public class MatrixOperationsTest extends AbstractTestCase {
	private static final Log log = LogFactory.getLog(CRFIOTest.class);
	boolean debug = log.isDebugEnabled();

	public void testMatrixExponential() throws Exception {
		//S = [-2.050871e-03, 5.984135e-04, 9.348390e-04, 5.176184e-04; 7.120564e-04, -2.164514e-03, 9.348390e-04, 5.176184e-04; 7.120564e-04, 5.984135e-04, -1.828088e-03, 5.176184e-04; 7.120564e-04, 5.984135e-04, 9.348390e-04, -2.245309e-03]

			System.out.println("Hello world");

		
		double[][] X = new double[4][4];
		X[0][0]=-2.050871e-03;   X[0][1]= 5.984135e-04;   X[0][2]= 9.348390e-04;   X[0][3]= 5.176184e-04;
		X[1][0]= 7.120564e-04;   X[1][1]=-2.164514e-03;   X[1][2]= 9.348390e-04;   X[1][3]= 5.176184e-04;
		X[2][0]= 7.120564e-04;   X[2][1]= 5.984135e-04;   X[2][2]=-1.828088e-03;   X[2][3]= 5.176184e-04;
		X[3][0]= 7.120564e-04;   X[3][1]= 5.984135e-04;   X[3][2]= 9.348390e-04;   X[3][3]=-2.245309e-03;

		DoubleMatrix2D S = new DenseDoubleMatrix2D(4,4);
		S.assign(X);
		System.out.println("The matrix S is: ");
		System.out.println(ColtUtil.format(S));
		
		DoubleMatrix2D T = new DenseDoubleMatrix2D(4,4);
		ColtUtil.exponentiate_real_matrix(S,T,200);
		System.out.println("The matrix T (the computed exponential of S) is: ");
		System.out.println(ColtUtil.format(T));		
		
		
		// Below we paste in the correct answer as computed by Matlab U = expm(S)
		double[][] Y = new double[4][4];
		Y[0][0]=0.9980;  Y[0][1]=0.0006;  Y[0][2]=0.0009;  Y[0][3]=0.0005;
		Y[1][0]=0.0007;  Y[1][1]=0.9978;  Y[1][2]=0.0009;  Y[1][3]=0.0005;
		Y[2][0]=0.0007;  Y[2][1]=0.0006;  Y[2][2]=0.9982;  Y[2][3]=0.0005;
		Y[3][0]=0.0007;  Y[3][1]=0.0006;  Y[3][2]=0.0009;  Y[3][3]=0.9978;
		DoubleMatrix2D U = new DenseDoubleMatrix2D(4,4);		
		U.assign(Y);
		System.out.println("The matrix U (the computed exponential of S using Matlab's U = expm(S)) is: ");
		System.out.println(ColtUtil.format(U));			
		
		for (int i=0; i<4; i++) {
			for (int j=0; j<4; j++) {
				assertEquals(T.getQuick(i,j),U.getQuick(i,j),0.001);
			}
		}
	    
	}
	

}
