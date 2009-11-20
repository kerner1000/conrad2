package calhoun.analysis.crf.executables.test;

import calhoun.analysis.crf.Conrad;


public class DaveTest {

	public static void main(String[] args) throws Exception {
		Conrad.main(new String[] { "train", "tricycle13.xml", "train_100_1", "baseline_cml_train_100_0_features.ser"});
		//Conrad.main(new String[] { "trainFeatures", "baseline_cml.xml", "train_100_1", "baseline_cml_train_100_0_features.ser"});
		//Conrad.main(new String[] { "trainWeights", "baseline_cml_train_100_0_features.ser", "train_100_1_sub", "baseline_cml_train_100_0.ser"});
		//Node marginals at seq 10 last position (1517): 5.203959e-114, 0.000000e+00, 0.000000e+00, 1.172256e-41, 1.000000e+00, 1.471485e-22, 8.193559e-10, 0.000000e+00, 0.000000e+00, 0.000000e+00, 2.835172e-47, 4.451058e-133, 7.265315e-96
		//Conrad.main(new String[] { "train", "z:/CRF/.xml", "sample/fly", "test/working/flyModel.ser"});
		//Conrad.main(new String[] { "train", "sample/fly/modifiedModelAllLengths.xml", "sample/fly", "test/working/flyModel.ser"});
		//Conrad.main(new String[] { "test", "test/working/flyModel.ser", "sample/fly", "test/working/flyOut"});
		//CRFTester.main(new String[] {"test/input/cryptoGeneSet", "test/working/cryptoGeneSet"});
		//CRFTester.main(new String[] {"test/input/crypto5wayGeneSet", "test/working/crypto5wayGeneSet"});
		//CRFTester.main(new String[] {"test/input/cryptoAOF", "test/working/cryptoAOF"});
		//CRFTester.main(new String[] {"test/input/cryptoAOFDebug", "test/working/cryptoAOFDebug"});
	}
}
