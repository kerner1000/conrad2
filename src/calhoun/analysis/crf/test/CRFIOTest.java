package calhoun.analysis.crf.test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.Conrad;
import calhoun.analysis.crf.executables.InputSequenceSubsetter;
import calhoun.analysis.crf.executables.PartitionRandomlyTrainTestFiles;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.util.AbstractTestCase;
import calhoun.util.FileUtil;

public class CRFIOTest extends AbstractTestCase {
	private static final Log log = LogFactory.getLog(CRFIOTest.class);
	boolean debug = log.isDebugEnabled();
	
	
	public void testArrayReaderInt() throws Exception {

		String fileName = "test/input/test_tabbed_array_reader_int.txt";		
		String[][] s = FileUtil.readFlatFile(fileName);
		
		assertEquals(s[0][0],"14");
		assertEquals(s[0][1],"276470");			
	}

	public void testInputSequenceReader() throws Exception {
		String configFile = "test/input/test_subsetting_configfile.xml";
		String inputFile = "test/input/test_subsetting_inputfile.txt";
		String outputFile = "test/working/test_subsetting_output.txt";
		
		// Define a model manager from file
		Conrad c = new Conrad(configFile);
		
		// Read a list of inputsequences from file, of the type expected by model manager cm
		List<? extends TrainingSequence<?>> t = c.getInputHandler().readTrainingData(inputFile);
		assertEquals(t.size(),1);
		assertEquals(t.get(0).length(),7560);
		
		// Take subsets of the first of the InputSequences you read
		List<TrainingSequence<?>> s = new ArrayList<TrainingSequence<?>>();
		s.add(t.get(0).subSequence(1,10));
		s.add(t.get(0).subSequence(1,20));
		
		// Write these subsetted InputSequences to file; they should be of exactly the same type as before.
		c.getInputHandler().writeTrainingData(outputFile, s);

		// NOT AUTOMATABLE: Inspect the file by hand make sure it is what you think.
		
		// Can the file we just wrote be read in again?
		List<? extends TrainingSequence<?>> u = c.getInputHandler().readTrainingData(outputFile);
		assertEquals(2, u.size());
		assertEquals(10, u.get(0).length());
		assertEquals(20, u.get(1).length());
	}
	
	public void testInputSequenceFileReader() throws Exception {
		String configFile = "test/input/inputFilesTest/test_subsetting_configfile.xml";
		String configFile2 = "test/input/inputFilesTest/test_subsetting_configfile2.xml";
		String inputFile = "test/input/inputFilesTest";
		String outputFile = "test/working/test_file_subsetting_output.txt";
		String matchFile = "test/output/test_file_subsetting_output.txt";
		
		// Define a model manager from file
		Conrad c = new Conrad(configFile);
		
		// Read a list of inputsequences from file, of the type expected by model manager cm
		List<? extends TrainingSequence<?>> t = c.getInputHandler().readTrainingData(inputFile);
		assertEquals(2, t.size());
		assertEquals(7560, t.get(0).length());
		
		// Take subsets of the first of the InputSequences you read
		List<TrainingSequence<?>> s = new ArrayList<TrainingSequence<?>>();
		s.add(t.get(0).subSequence(1,10));
		s.add(t.get(0).subSequence(1,20));
		s.add(t.get(1).subSequence(1,10));
		
		// Write these subsetted InputSequences to file; they should be of exactly the same type as before.
		Conrad c2 = new Conrad(configFile2);
		c2.getInputHandler().writeTrainingData(outputFile, s);

		assertFilesMatch(matchFile, outputFile);
	}
	
	
	public void testInputSequenceSubsetterCommandLine() throws Exception {
		String configFile = "test/input/test_subsetting_configfile.xml";
		String inputFile = "test/input/test_subsetting_inputfile.txt";
		String regionsFile = "test/input/test_subsetting_regions.txt";
		String outputFile = "test/working/test_subsetting_output.txt";
		
		// Use a command line program to take specified subset regions
		InputSequenceSubsetter.main(new String[] {configFile,inputFile,regionsFile,"5",outputFile,"0"});
		
		// make sure you get what you expected
		Conrad c = new Conrad(configFile);
		List<? extends TrainingSequence<?>> u = c.getInputHandler().readTrainingData(outputFile);
		assertEquals(u.size(),2);
		assertEquals(u.get(0).length(),15);
		assertEquals(u.get(1).length(),30);
	}
	
	public void testInputSequenceSubsetterCommandLineSeparateFiles() throws Exception {
		String configFile = "test/input/inputFilesTest/test_subsetting_configfile.xml";
		String inputFile = "test/input/inputFilesTest";
		String regionsFile = "test/input/test_subsetting_regions.txt";
		String outputFile = "test/working/test_subsetting_output_files";
		
		File f = new File(outputFile);
		f.mkdirs();
		
		// Use a command line program to take specified subset regions
		InputSequenceSubsetter.main(new String[] {configFile,inputFile,regionsFile,"5",outputFile,"0"});
		
		// make sure you get what you expected
		assertFilesMatch("test/working/test_subsetting_output_files/hidden.dat", "test/output/test_subsetting_output_files/hidden.dat");
		assertFilesMatch("test/working/test_subsetting_output_files/ref.dat", "test/output/test_subsetting_output_files/ref.dat");
		assertFilesMatch("test/working/test_subsetting_output_files/name.dat", "test/output/test_subsetting_output_files/name.dat");
		assertFilesMatch("test/working/test_subsetting_output_files/aln.dat", "test/output/test_subsetting_output_files/aln.dat");
	}
	
	public void testInputSequenceSubsetterCommandLineForcingGenic() throws Exception {
		String configFile = "test/input/test_subsetting_configfile.xml";
		String inputFile = "test/input/test_subsetting_inputfile.txt";
		String regionsFile = "test/input/test_subsetting_regions.txt";
		String outputFile = "test/working/test_subsetting_output.txt";
		
		// Use a command line program to take specified subset regions
		InputSequenceSubsetter.main(new String[] {configFile,inputFile,regionsFile,"5",outputFile,"1"});
		
		// make sure you get what you expected
		Conrad c = new Conrad(configFile);
		List<? extends TrainingSequence<?>> u = c.getInputHandler().readTrainingData(outputFile);
		assertEquals(u.size(),0);
	}
	
	
	public void testTrainTestSplitCommandLine() throws Exception {
		String configFile = "test/input/test_subsetting_configfile.xml";
		String inputFile = "test/input/test_traintest_split.txt";
		String nTrainArg = "1";
		String outputTrain = "test/working/test_subsetting_output.txt";
		String outputTest = "test/working/test_subsetting_output.txt";		
		
		// Use a command line program to take specified subset regions
		PartitionRandomlyTrainTestFiles.main(new String[] {configFile,inputFile,nTrainArg,outputTrain,outputTest});
		
		// make sure you get what you expected
		Conrad c = new Conrad(configFile);
		List<? extends TrainingSequence<?>> u = c.getInputHandler().readTrainingData(outputTest);
		assertEquals(u.size(),2);
	}
	
}
