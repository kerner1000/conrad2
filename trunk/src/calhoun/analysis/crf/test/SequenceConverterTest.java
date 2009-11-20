package calhoun.analysis.crf.test;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.executables.HiddenSequenceConverter;
import calhoun.analysis.crf.io.InputHandler;
import calhoun.analysis.crf.io.InputHandlerInterleaved;
import calhoun.analysis.crf.io.SequenceConverter;
import calhoun.analysis.crf.io.StringInput;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.util.AbstractTestCase;
import calhoun.util.Assert;

public class SequenceConverterTest  extends AbstractTestCase {
	
	private static final Log log = LogFactory.getLog(SequenceConverterTest.class);
	boolean debug = log.isDebugEnabled();

	
	public void testTricycle13Interval13Conversion() throws Exception {
		System.out.println("Testing Hidden Sequence Converter between Interval13 and Tricycle13 and back.");
		
		//		 This has nothing to do with GTF but this config file does what we want
		//ConfigurableModelManager cm= new ConfigurableModelManager("test/input/configGTF.txt");  

		InputHandler ih = new InputHandlerInterleaved(new StringInput());
		List<? extends TrainingSequence<?>> data = ih.readTrainingData("test/input/shortTrain.tricycle13.txt");
		List<TrainingSequence<Character>> listOfSeqs = new ArrayList<TrainingSequence<Character>>();		
		
		convertSequence(data, listOfSeqs);
		
		ih.writeTrainingData("test/working/shortTrain.interval13.txt", listOfSeqs);
		
		// Test below is a little weak right now (20060705) as I used this test to create the file in input; but will get
		// stronger as I plan to use the file "test/input/shortTrain.interval13.txt" elsewhere in testing 
		assertFilesMatch("test/input/interval13/data/shortTrain.interval13.txt","test/working/shortTrain.interval13.txt");
		
		//
		// SAME TEST, DIFFERENT FILE
		//
		ih = new InputHandlerInterleaved(new StringInput());
		List<? extends TrainingSequence<?>> data2 = ih.readTrainingData("test/input/negStrandGene.tricycle13.txt");
		List<TrainingSequence<Character>> listOfSeqs2 = new ArrayList<TrainingSequence<Character>>();		
		
		convertSequence(data2, listOfSeqs2);
		
		ih.writeTrainingData("test/working/negStrandGene.interval13.txt", listOfSeqs2);
		
		assertFilesMatch("test/input/interval13/data/negStrandGene.interval13.txt","test/working/negStrandGene.interval13.txt");
	}
	
	private void convertSequence(List<? extends TrainingSequence<?>> data, List<TrainingSequence<Character>> listOfSeqs) {
		for (int i=0; i<data.size(); i++) {
			TrainingSequence<Character> seq = (TrainingSequence<Character>) data.get(i);
			int len = seq.length();
			
			int[] oldy = new int[len];
			for (int j=0; j<len; j++) { oldy[j] = seq.getY(j); }
		
			SequenceConverter.convertSeqFromTricycle13ToInterval13(seq);
			
			listOfSeqs.add(seq);
			
			SequenceConverter.convertSeqFromInterval13ToTricycle13(seq);			

			int[] newy = new int[len];
			for (int j=0; j<len; j++) {
				newy[j] = seq.getY(j);
				Assert.a(newy[j] == oldy[j]);
			}
			
			SequenceConverter.convertSeqFromTricycle13ToInterval13(seq);
		}
	}
	
	
	public void testSequenceConverter() throws Exception {
		System.out.println("Testing Hidden Sequence Converter.");
		String[] args13 = new String[4];
		String[] args39 = new String[4];
		
		args13[0] = "13to39";
		args13[1] = "test/input/ghmmUsingCrfIO/basic_ModelDSComp.xml"; 	// config file
		args13[2]   = "test/input/ghmmUsingCrfIO/cnDT_chr14_60k.txt";	
		args13[3] = "test/working/cnDT_chr14_60k_39states.txt";

		args39[0] = "39to13";
		args39[1] = "test/input/ghmmUsingCrfIO/basic_ModelDSComp.xml"; 	// config file
		args39[2] = args13[3];
		args39[3] = "test/working/cnDT_chr14_60k_13states.txt";
		
		HiddenSequenceConverter.main(args13);
		HiddenSequenceConverter.main(args39);
	}

	public void testSequenceConverter13models() throws Exception {
		System.out.println("Testing Hidden Sequence Converter.");
		String[] argsTri2Int = new String[3];
		String[] argsInt2Tri = new String[3];
		
		argsTri2Int[0] = "tricycle13ToInterval13";
		argsTri2Int[1] = "test/input/interval13/data/splitInputOneGeneTrain/hidden.tricycle13.dat";	
		argsTri2Int[2] = "test/working/shortTrain.interval13.txt";

		argsInt2Tri[0] = "interval13ToTricycle13";
		argsInt2Tri[1] = "test/input/interval13/data/splitInputOneGeneTrain/hidden.dat";	
		argsInt2Tri[2] = "test/working/shortTrain.tricycle13.txt";
		
		HiddenSequenceConverter.main(argsTri2Int);
		HiddenSequenceConverter.main(argsInt2Tri);
		
		assertFilesMatch(argsTri2Int[1],argsInt2Tri[2]);
		assertFilesMatch(argsTri2Int[2],argsInt2Tri[1]);
	}
	
	
	
}
