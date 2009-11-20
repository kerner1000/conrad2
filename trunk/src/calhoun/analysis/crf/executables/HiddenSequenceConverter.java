package calhoun.analysis.crf.executables;


import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import calhoun.analysis.crf.Conrad;
import calhoun.analysis.crf.io.SequenceConverter;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.util.Assert;
import calhoun.util.FileUtil;


/** Converts hidden sequences between varioud models.
 * <p>
 * Example:<p> 
 * <code>HiddenSequenceConverter tricycle13ToInterval13 <config file> <input file> <output file></code>
 * 
 * <b>HiddenSequenceConverter 13to39</b>converts a 13 state sequence to a 39 state sequence, and 
 * produces a file that contains a hidden sequence of values [0, 38].
 * <b>HiddenSequenceConverter 39to13</b> converts a 39 state sequence to a 13 state sequence,
 * and produces a file that contains a hidden sequence of values [0, 12].
 */
public class HiddenSequenceConverter {
	
	private static void explainUsage() {
		System.err.println("Usage error.  Call must be the following format:");
		System.err.println("> HiddenSequenceConverter 13to39 <config file> <input file> <output file>");
		System.err.println("> HiddenSequenceConverter 39to13 <config file> <input file> <output file>");
		System.err.println("> HiddenSequenceConverter tricycle13ToInterval13 <config file> <input file> <output file>");
		System.err.println("> HiddenSequenceConverter interval13ToTricycle13 <config file> <input file> <output file>");
		System.err.println("> HiddenSequenceConverter tricycle13ToInterval13 <input file.dat> <output file.dat>");
		System.err.println("> HiddenSequenceConverter interval13ToTricycle13 <input file.dat> <output file.dat>");
		
		System.exit(0);		
	}
	
	public static void main(String[] args) throws IOException 
	{		
		if ((args.length < 3) || (args.length>4)) {
			System.out.println("ERROR - number of args is less than 3 or greater than 4");
			explainUsage(); }
		
		if (args.length==3) {
			
			System.out.println("Hidden Sequence Converter");
			
			String convertType = args[0];
			String inputFile   = args[1];
			String outputFile  = args[2];
			
			String[] data = FileUtil.readLines(inputFile);
			
			BufferedWriter  w = new BufferedWriter(new FileWriter(outputFile));
			
			System.out.println("Size of data is " + data.length );
			
			for (int j=0; j<data.length; j++) {
				if (convertType.equals("tricycle13ToInterval13")) {
					w.write( SequenceConverter.convertSeqFromTricycle13ToInterval13(data[j]) + "\n");
				}
				else if (convertType.equals("interval13ToTricycle13")) {
					w.write( SequenceConverter.convertSeqFromInterval13ToTricycle13(data[j]) + "\n");
				}
				else {
					System.out.println("ERROR - convertType is not recognized");
					explainUsage();
					Assert.a(false,"convertType = " + convertType);
				}	
			}
			w.close();
		}
		
		if (args.length==4) {
			System.out.println("Hidden Sequence Converter");
			
			String convertType = args[0];
			String configFile  = args[1];
			String inputFile   = args[2];
			String outputFile  = args[3];
			int i;
			
			// Define a model manager from file
			Conrad crf = new Conrad(configFile);
			
			// Read a list of input sequences from file, of the type expected by model manager cm
			List<? extends TrainingSequence<?>> data = crf.getInputHandler().readTrainingData(inputFile);
			List<TrainingSequence<Character>> listOfSeqs = new ArrayList<TrainingSequence<Character>>();		
			
			
			for (i=0; i<data.size(); i++)
			{
				TrainingSequence<Character> seq = (TrainingSequence<Character>) data.get(i);
				System.out.println("Converting sequence " + i + ".");
				if (convertType.equals("13to39")) {
					SequenceConverter.convertSeqFrom13To39(seq.getTrainingComponent("ref"));
				}
				else if (convertType.equals("39to13")) {
					SequenceConverter.convertSeqFrom39To13(seq.getTrainingComponent("ref"));
				}
				else if (convertType.equals("tricycle13ToInterval13")) {
					SequenceConverter.convertSeqFromTricycle13ToInterval13(seq.getTrainingComponent("ref"));
				}
				else if (convertType.equals("interval13ToTricycle13")) {
					SequenceConverter.convertSeqFromInterval13ToTricycle13(seq.getTrainingComponent("ref"));
				}
				else {
					explainUsage();
					Assert.a(false,"convertType = " + convertType);
				}
				listOfSeqs.add(seq);
			}
			
			// Output new file.
			crf.getInputHandler().writeTrainingData(outputFile, listOfSeqs);
		}
		
		System.out.println("DONE");
	}
	
}
