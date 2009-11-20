package calhoun.analysis.crf.executables;

import java.io.File;
import java.io.FileWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang.StringUtils;

import calhoun.analysis.crf.Conrad;
import calhoun.util.FileUtil;

/** Tests a CRF using a variety of train and test sets and a variety of different models.
 * Takes an input directory and an output directory on the command line.
 * <pre> 
 * Makes an output directory 
 * xxxModel yyyTrain zzzTest 
 * 
 * Saves a model file xxxyyyModel.ser
 * All weights go in a weight table - xxxyyy w1 w2 w3 w4 ...
 * Puts all results in a result table
 * xxxyyy #features trainTime trainPerformance test1 test2 test3 test4
 * </pre> 
 */
public class CRFTester {
	public static void main(String[]args) throws Exception {
		String dir = args[0];
		String outputDir = args[1];
		List<String> modelFiles = new ArrayList();
		List<String> trainingSets = new ArrayList();
		List<String> testSets = new ArrayList();
		File f = new File(dir);
		File[] names = f.listFiles();
		for(File name : names) {
			String fName = name.getAbsolutePath();
			if(fName.indexOf("Model") != -1) {
				modelFiles.add(fName);
			}
			else if(fName.indexOf("Train") != -1) {
				trainingSets.add(fName);
			}
			else if(fName.indexOf("Test") != -1) {
				testSets.add(fName);
			}
		}
		runTests(modelFiles, trainingSets, testSets, outputDir);
	}

	static void runTests(List<String> modelFiles, List<String> trainingSets, List<String> testSets, String outputDir) throws Exception {
		File dir = new File(outputDir);
		File results = new File(dir, "results.csv");
		File weights = new File(dir, "weights.csv");
		dir.mkdirs();
		Writer rw = new FileWriter(results);
		rw.write("Model\tTraining\tFeatures\tTrainTime\tTrain\t"+StringUtils.join(testSets.iterator(),'\t'));
		Writer ww = new FileWriter(weights);
		for(String model : modelFiles) {
			String modelName = FileUtil.getBaseAndExtension(new File(model))[0];
			for(String train : trainingSets) {
				String trainName = FileUtil.getBaseAndExtension(new File(train))[0];
				//Conrad r = new Conrad(model);
				Conrad r = new Conrad(model);
				System.out.println("A model was successfully initiailized");
				// Train
				r.train(train);
				System.out.println("Model has been trained");
				// Write model
				r.write(new File(dir, modelName+trainName+".ser").getAbsolutePath());
				System.out.println("Trained model has been written to file");
				// Write weights
				ww.write(modelName+trainName+'\n'+r.printWeights());
				System.out.println("Weights written to file");
				
				// Write results
				System.out.println("Testing performance of model on training data");
				r.test(train, new File(dir, "predicted"+modelName+trainName+trainName+".txt").getAbsolutePath());
				//System.out.println(modelName+'\t'+trainName+'\t'+r.getModel().getNumFeatures()+'\t'+String.format("%.2f",ret.getAccuracy()*100));
				//rw.write(modelName+'\t'+trainName+'\t'+r.getModel().getNumFeatures()+'\t'+String.format("%.2f",ret.getAccuracy()*100));
				
				for(String test : testSets) {
					String testName = FileUtil.getBaseAndExtension(new File(test))[0];
					r.test(test, new File(dir, "predicted"+modelName+trainName+testName+".txt").getAbsolutePath());
					//rw.write('\t'+String.format("%.2f",ret.getAccuracy()*100));
					//ret.writeData();
				}
				rw.append("\n");
				System.out.println("Model has been tested on all data; moving to next model or stopping");
			}
		}
		rw.close();
		ww.close();
		System.out.println("Done with Run Tests");
	}
}
