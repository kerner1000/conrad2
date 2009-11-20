package calhoun.analysis.crf.executables;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import calhoun.analysis.crf.Conrad;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.util.Assert;

/** takes a set of sequences and randomly divides them into training adn test sets.
 */
public class PartitionRandomlyTrainTestFiles {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws Exception {
		String configFile  = args[0];  // "test/input/test_subsetting_configfile.xml";
		String inputFile   = args[1];  // "test/input/test_subsetting_inputfile.txt";  
		int    nTrain      = Integer.parseInt(args[2]);
		String outputTrain = args[3];  // "test/output/test_subsetting_output.txt";
		String outputTest  = args[4];  // "test/output/test_subsetting_output.txt";
			
		// Define a model manager from file
		Conrad c = new Conrad(configFile);
		
		// Read a list of inputsequences from file, of the type expected by model manager cm
		Iterator<? extends TrainingSequence<?>> iter = c.getInputHandler().readTrainingData(inputFile).iterator();
		List<TrainingSequence> t = new ArrayList<TrainingSequence>();
		while(iter.hasNext()) {
			t.add(iter.next());
		}
		
		int nTotal = t.size();
		Assert.a(nTrain > 0);
		Assert.a(nTrain < nTotal);
		boolean[] trainflag = randomBooleanVector(nTotal,nTrain);

		
		List<TrainingSequence<?>> train = new ArrayList<TrainingSequence<?>>();
		List<TrainingSequence<?>> test  = new ArrayList<TrainingSequence<?>>();
		
		for (int i=0; i<nTotal; i++) {
			if (trainflag[i]) {
				train.add(t.get(i));
			} else {
				test.add(t.get(i));
			}
		}

		c.getInputHandler().writeTrainingData(outputTrain, train);
		c.getInputHandler().writeTrainingData(outputTest, test);
	}

	private static boolean[] randomBooleanVector(int total, int train) {
		boolean[] bv = new boolean[total];
		for (int i=0; i<total; i++) { bv[i] = false; }		
		for (int i=0; i<train; i++) {
			Random r = new Random();
			int x;
			do {
				x = r.nextInt(total);
			} while (bv[x]);
			bv[x] = true;
		}
		
		return bv;
	}
}
