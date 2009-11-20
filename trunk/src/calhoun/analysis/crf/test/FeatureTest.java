package calhoun.analysis.crf.test;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.List;

import org.apache.commons.lang.StringUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.Conrad;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.features.tricycle13.KmerFeatures;
import calhoun.analysis.crf.io.InputHandlerInterleaved;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.analysis.crf.io.StringInput;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.analysis.crf.solver.check.ArrayFeatureList;
import calhoun.util.AbstractTestCase;
import calhoun.util.Assert;

public class FeatureTest extends AbstractTestCase {
	private static final Log log = LogFactory.getLog(FeatureTest.class);
	boolean debug = log.isDebugEnabled();
	
	public void testKmerFeatures() throws Exception {
		List<TrainingSequence<Character>> data = (List<TrainingSequence<Character>>) new InputHandlerInterleaved(new StringInput()).readTrainingData("test/input/zeroOrderTest.txt");
		KmerFeatures kf = new KmerFeatures();
		kf.setRareThreshold(0);
		kf.train(0, new ZeroOrderManager(), data);
		assertEquals(2, kf.getNumFeatures());
		assertEquals(expectedProb(53, 30), kf.getKmerProb(0, "A", 0), .01);
		assertEquals(expectedProb(30, 53), kf.getKmerProb(0, "A", 1), .01);
		assertEquals(expectedProb(55, 38), kf.getKmerProb(0, "T", 0), .01);
		assertEquals(expectedProb(38, 55), kf.getKmerProb(0, "T", 1), .01);
		assertEquals(expectedProb(21, 91), kf.getKmerProb(0, "C", 0), .01);
		assertEquals(expectedProb(91, 21), kf.getKmerProb(0, "C", 1), .01);
		assertEquals(expectedProb(23, 82), kf.getKmerProb(0, "G", 0), .01);
		assertEquals(expectedProb(82, 23), kf.getKmerProb(0, "G", 1), .01);
	}

	double expectedProb(int a, int b) {
		return Math.log((a+1)/(float)(a+b+2));
	}
	
	public void testFeatureNames() throws Exception {
		List<? extends TrainingSequence<?>> data = new InputHandlerInterleaved(new StringInput()).readTrainingData("test/input/zeroOrderTest.txt");
		Conrad runner = ZeroOrderManager.getCRF();
		runner.trainFeatures(data);
		assertEquals("Start.lowGC", runner.getModel().getFeatureName(0));
		assertEquals("End.highGC", runner.getModel().getFeatureName(3));
		assertEquals("Kmer.lowGC.0", runner.getModel().getFeatureName(4));
		assertEquals("Kmer.highGC.0", runner.getModel().getFeatureName(5));
		assertEquals("Edge.lowGC-lowGC", runner.getModel().getFeatureName(6));
		assertEquals(10, runner.getModel().getNumFeatures());
	}

	public void testWriteFeatures() throws Exception {
		List<? extends TrainingSequence<?>> data = new InputHandlerInterleaved(new StringInput()).readTrainingData("test/input/zeroOrderTest.txt");
		Conrad runner = ZeroOrderManager.getCRF();
		runner.trainFeatures(data);
		writeFeatures(runner, "test/working/zeroOrderFeatures.csv", data.get(0));
		assertFilesMatch("test/output/zeroOrderFeatures.csv", "test/working/zeroOrderFeatures.csv");
	}

	public void testComponentFeatures() throws Exception {
		// Test that component features work correctly
		Conrad runner = new Conrad("test/input/componentFeatures.xml");
		List<? extends TrainingSequence<?>> data = runner.getInputHandler().readTrainingData("test/input/testTrain.txt");
		runner.trainFeatures(data);
		ArrayFeatureList result = new ArrayFeatureList(runner.getModel()); 
		result.evaluateNode(data.get(0), 0, 0);
	}

	public void writeFeatures(Conrad runner, String file, InputSequence data) throws IOException {
		ModelManager model = runner.getModel();
		boolean training = TrainingSequence.class.isInstance(data);
		int numFixed = 2 + (training ? 1 : 0);
		int totalFeatures = model.getNumFeatures();
		String[] line = new String[numFixed + totalFeatures];
		
		Writer fw = new BufferedWriter(new FileWriter(file));
		// Write the header line
		line[0] = "pos";
		line[1] = "x";
		if(training) {
			line[2] = "y";
		}
		for (int i = 0; i < totalFeatures; ++i) {
			line[i + numFixed] = model.getFeatureName(i);
		}
		fw.write('#' + StringUtils.join(line, "\t") + '\n');
		
		for (int i = 0; i < data.length(); ++i) {
			line[0] = Integer.toString(i);
			line[1] = data.getX(i).toString();
			if(training) {
				line[2] = Integer.toString(((TrainingSequence) data).getY(i));
			}
			for(int j = 0; j<totalFeatures; ++j) {
				line[numFixed+j] = "0";
			}
			int nStates = model.getNumStates();
			ArrayFeatureList result = new ArrayFeatureList(model);
			for(int state = 0; state < nStates; ++state) {
				result.evaluateNode(data, i, state);
				updateFeatureValues(result, line, numFixed);
				if(i > 0) {
					for(int prevState = 0; prevState < nStates; ++prevState) {
						result.evaluateEdge(data, i, prevState, state);
						updateFeatureValues(result, line, numFixed);
					}
				}
			}
			fw.write(StringUtils.join(line, "\t") + '\n');
		}
		fw.close();
	}

	void updateFeatureValues(ArrayFeatureList result, String[] line, int numFixed) {
		int[] indices = result.getIndices();
		double[] vals = result.getValues();
		int size = result.size();
		for(int i = 0; i<size; ++i) {
			Assert.a(line[numFixed + indices[i]].equals("0"), "Feature had a previous value: "+line[numFixed + indices[i]]);
			line[numFixed + indices[i]] = Float.toString((float)vals[i]);
		}
	}
}
