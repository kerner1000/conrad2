package calhoun.analysis.crf.executables.test;

import calhoun.analysis.crf.Conrad;
import calhoun.util.DenseBooleanMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

/** Not really a test case, just put here because it is a convenient function to run the profiler on to check speed */
public class ProfilerRun {

	public static void main(String[] args) throws Exception {
		semiMarkovTraining();
		//semiMarkovViterbi();
		
		//testArrayAccess();
	}
	
	static void semiMarkovTraining() throws Exception {
		Conrad.main(new String[] {"train", "sample/baseline.xml", "train_100_1", "test/working/profilerModel.ser"});
	}
	
	static void semiMarkovViterbi() throws Exception {
		Conrad.main(new String[] {"test", "test/working/interval13BaselineModelTest.ser", "test/input/interval13/data/oneGeneTruncatedTrain.interval13.txt", "test/working/interval13SemiMarkovModelTestPredicted.txt"} );
	}
	
	static void testArrayAccess() {
		System.out.println("Allocating");
		boolean[] array= new boolean[500000];
		DenseBooleanMatrix2D bArray= new DenseBooleanMatrix2D(5, 100000);

		double[] dArray= new double[100000];
		DenseDoubleMatrix2D dArray2= new DenseDoubleMatrix2D(10, 10000);
		
		
		System.out.println("Begin");
		for(int i=0; i<1000; ++i) {
			test1(array);
			test1b(bArray);
			test2(dArray);
			test2b(dArray2);
		}
	}

	static void test1(boolean[] array) {
		int count = 0;
		for(int i = 0; i<5; i++) {
			for(int j = 0; j<100000; j++) {
				count += array[i*100000+j] ? 1 : 0;
			}
		}
	}
	static void test1b(DenseBooleanMatrix2D array) {
		int count = 0;
		for(int i = 0; i<5; i++) {
			for(int j = 0; j<100000; j++) {
				count += array.getQuick(i, j) ? 1 : 0;
			}
		}
	}
	static void test2(double[] array) {
		int count = 0;
		for(int i = 0; i<10; i++) {
			for(int j = 0; j<10000; j++) {
				count += array[j*10 + i]==0 ? 1 : 0;
			}
		}
	}
	static void test2b(DenseDoubleMatrix2D array) {
		int count = 0;
		for(int i = 0; i<10; i++) {
			for(int j = 0; j<10000; j++) {
				count += array.getQuick(i, j)==0 ? 1 : 0;
			}
		}
	}
}
