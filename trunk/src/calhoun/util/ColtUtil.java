package calhoun.util;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;

import cern.colt.function.DoubleFunction;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.EigenvalueDecomposition;

public class ColtUtil {

	public static void exponentiate_real_matrix(DoubleMatrix2D X, DoubleMatrix2D Y, int degree)
	// Sets Y to equal exp(X); X is assumed to be real.
	// Y can be expressed as a convergent power series in X:
	// Y = I + X + X*X/2 + X*X*X/6 + X*X*X*X/24 + ...
	// the power series is truncated after the X^(degree) term.
	// Since X is real and each term above is real, Y is also real.
	// An eigenvector of X with eigenvalue lambda is an eigenvector of Y with eigenvalue exp(lambda)
	{
		// S is a reall symmetric matrix, so it has real eigenvalues.
		/*
		EigenvalueDecomposition E = new EigenvalueDecomposition(X);
		
		Assert.a(degree >= 1);
		
		// A = V*D*V'
		DoubleMatrix2D V = E.getV();
		DoubleMatrix2D D = E.getD();			
		
		System.out.println("The matrix D is: ");
		System.out.println(ColtUtil.format(D));
		System.out.println("The matrix V is: ");
		System.out.println(ColtUtil.format(V));

		System.out.println("Real eigenvalues D is: ");
		System.out.println(ColtUtil.format(E.getRealEigenvalues().toArray()));
		System.out.println("Imaginary eigenvalues D is: ");
		System.out.println(ColtUtil.format(E.getImagEigenvalues().toArray()));
		*/
		/*
		for (int i=0; i<4; i++) {
			// only the diagonal entries:
			D.setQuick(i,i,Math.exp(D.getQuick(i,i)));
		}
		*/
		/*
		System.out.println("After exponentiating diagonal, The matrix D is: ");
		System.out.println(ColtUtil.format(D));
		
		DoubleMatrix2D B = new DenseDoubleMatrix2D(4,4); B.assign(0);
		V.zMult(D,B); // B=V*D;
		
		B.zMult(V,Y,1.0,0.0,false,true); // Y = 1.0*B*V' + 0.0*Y;
		*/
		
		//ColtMatrix Z;
		
		DenseDoubleMatrix2D EYE = new DenseDoubleMatrix2D(4,4);
		EYE.assign(0.0);  for (int i=0; i<4; i++) { EYE.setQuick(i,i,1.0); }  // Intially sets Y to the identity matrix; 
		Y.assign(EYE); 
		DenseDoubleMatrix2D Z = new DenseDoubleMatrix2D(4,4);   Z.assign(EYE);
		DenseDoubleMatrix2D W = new DenseDoubleMatrix2D(4,4);
		for (int j=1; j<=degree; j++) {
			W.assign(Z);
			W.zMult(X,Z,1.0/((double) j),0.0,false,false);   // Z = (1/j)*W*X + 0*Z = X^j/j!
			Z.zMult(EYE,Y,1.0,1.0,false,true);   // Y = Z*I + Y = Z+Y = Y+X^j/j!
		}
		
		//Y = ((ColtMatrix) X).exponential(200);
	
		/*public DoubleMatrix2D zMult(DoubleMatrix2D B,
                            DoubleMatrix2D C,
                            double alpha,
                            double beta,
                            boolean transposeA,
                            boolean transposeB)
Linear algebraic matrix-matrix multiplication; C = alpha * A x B + beta*C. C[i,j] = alpha*Sum(A[i,k] * B[k,j]) + beta*C[i,j], k=0..n-1. 
Matrix shapes: A(m x n), B(n x p), C(m x p). 
Note: Matrix shape conformance is checked after potential transpositions.
*/
		
		return;
	}
	
	
	public static void exponentiate_symmetric_matrix(DoubleMatrix2D X, DoubleMatrix2D Y)
	// Sets Y to equal exp(X); X is assumed to be real and symmetric.
	{
		// S is a reall symmetric matrix, so it has real eigenvalues.
		EigenvalueDecomposition E = new EigenvalueDecomposition(X); 
		
		// A = V*D*V'
		DoubleMatrix2D V = E.getV();
		DoubleMatrix2D D = E.getD();			
		
		for (int i=0; i<4; i++) {
			// only the diagonal entries:
			D.setQuick(i,i,Math.exp(D.getQuick(i,i)));
		}
		
		DoubleMatrix2D B = new DenseDoubleMatrix2D(4,4); B.assign(0);
		V.zMult(D,B); // B=V*D;
		
		B.zMult(V,Y,1.0,0.0,false,true); // Y = 1.0*B*V' + 0.0*Y;
	
		return;
	}
	
	public static double dotProduct(double[] a, double[] b) {
		int length = a.length;
		Assert.a(a.length == b.length, length, " length array dotted with array of length ", b.length);
		double ret = 0.0;
		for(int i=0; i<length; ++i) {
			ret += a[i]*b[i];
		}
		return ret;
	}
	
	public static String format(double[] str) {
		StringBuffer b = new StringBuffer();
		for (double d : str) {
			b.append(String.format("%e, ", d));
		}
		if(b.length() > 2) {
			b.setLength(b.length() - 2);
		}
		return b.toString();
	}

	public static void printToFile(double[] str, int numRows, int numCols) {
		try {
			Writer fout = new BufferedWriter(new FileWriter("test/working/crf_forwardPass_60kTEST.txt"));

			fout.write("Viterbi Map from CRF Forward Pass\n");
			fout.write("\n");
			fout.write("\t");
			for (int r=0; r<numRows; r++)
				fout.write(String.format("%1$11d", r));
			fout.write("\n\n");
			
			for (int c=0; c<numCols; c++)
			{
				fout.write(c + ")\t");
				for (int r=0; r<numRows; r++)
				{
					double val = str[numRows * c + r];
					if (Double.isInfinite(val))
						fout.write(String.format("%1$11.0f", 0.0f));
					else
						//fout.write(String.format("%e\t", bestScores[pos][st], 2));
						fout.write(String.format("%1$11.2f", val));
				}
				fout.write("\n");
			}
			fout.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}	
	}	
	
	public static void print(double[] str, int numRows, int numCols) {
		System.out.print("\n");
		for (int r=0; r<numRows; r++)
		{
			for (int c=0; c<numCols; c++)
			{
				double val = str[numRows * c + r];
				if (Double.isInfinite(val))
					val = 0.0;
				System.out.print(String.format("%e\t", val));
			}
			System.out.print("\n");
		}
	}
	
	public static void print(int[] str, int numRows, int numCols) {
		System.out.print("\n");
		for (int r=0; r<numRows; r++)
		{
			for (int c=0; c<numCols; c++)
			{
				System.out.print(String.format("%d\t", str[numRows * c + r]));
			}
			System.out.print("\n");
		}
	}
	
    public static double norm(double ar[]) {
        double v = 0;
        for (int f = 0; f < ar.length; f++)
            v += ar[f]*ar[f];
        return Math.sqrt(v);
    }

    public static String format(DoubleMatrix2D mi) {
		StringBuffer b = new StringBuffer();
		for (int i = 0; i < mi.rows(); ++i) {
			for (int j = 0; j < mi.columns(); ++j) {
				b.append(String.format("%e, ", mi.getQuick(i, j)));
			}
			b.setLength(b.length() - 2);
			b.append("; ");
		}
		return b.toString();
	}
	
	/** Returns the row with the largest value in the given column */
	public static int maxInColumn(DoubleMatrix2D m, int col) {
		double max = Double.NEGATIVE_INFINITY;
		int bestRow = -1;
		for(int j = 0; j<m.rows(); ++j) {
			double current = m.getQuick(j, col); 
			if(current > max) {
				max = current;
				bestRow = j;
			}
		}
		Assert.a(bestRow != -1);
		return bestRow;
	}

	public static int maxInColumn(double[] m, int numRows, int col) {
		double max = Double.NEGATIVE_INFINITY;
		int bestRow = -1;
		int offset = col*numRows;
		for(int j = 0; j<numRows; ++j) {
			double current = m[offset + j]; 
			if(current > max) {
				max = current;
				bestRow = j;
			}
		}
		Assert.a(bestRow != -1);
		return bestRow;
	}

	public static double maxValueInColumn(double[] m, int numRows, int col) {
		double max = Double.NEGATIVE_INFINITY;
		int offset = col*numRows;
		for(int j = 0; j<numRows; ++j) {
			double current = m[offset + j]; 
			if(current > max) {
				max = current;
			}
		}
		return max;
	}

	public static final DoubleFunction exp = new  DoubleFunction() {
		public double apply(double a) {
			return Math.exp(a);
		}
	};
	public static final DoubleFunction ln = new DoubleFunction() {
		public double apply(double a) {
			return Math.log(a);
		}
	};
	public static void scalarmultiply(DoubleMatrix2D S, final double t) {
		DoubleFunction mult = new DoubleFunction() {
			public double apply(double a) {
				return t*a;
			}
		};
		S.assign(mult);	
	}
	
}
