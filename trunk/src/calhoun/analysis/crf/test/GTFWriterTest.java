package calhoun.analysis.crf.test;

import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.Conrad;
import calhoun.analysis.crf.io.InputHandler;
import calhoun.analysis.crf.io.InputHandlerInterleaved;
import calhoun.analysis.crf.io.OutputHandlerGeneCallStats;
import calhoun.analysis.crf.io.StringInput;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.util.AbstractTestCase;

public class GTFWriterTest extends AbstractTestCase {
	private static final Log log = LogFactory.getLog(CacheTest.class);
	boolean debug = log.isDebugEnabled();
	
	public void testWriteGTF() throws Exception {
		Conrad c = new Conrad("test/input/configGTF.xml");

		// Test Tricycle13 model
		String gtfFile = "test/working/shortTrain.gtf";
		InputHandler ih = new InputHandlerInterleaved(new StringInput());
		List<? extends TrainingSequence<?>> data = ih.readTrainingData("test/input/shortTrain.tricycle13.txt");
		OutputHandlerGeneCallStats oh = new OutputHandlerGeneCallStats(c.getModel(), ih);
		oh.writeGTF(data, gtfFile);
		assertFilesMatch("test/output/shortTrain.gtf", gtfFile);
		
		gtfFile = "test/working/testGTF.gtf";
		List<? extends TrainingSequence<?>> data2 = c.getInputHandler().readTrainingData("test/input/testGTF.txt");
		oh = new OutputHandlerGeneCallStats(c.getModel(), c.getInputHandler());
		oh.writeGTF(data2, gtfFile);
		assertFilesMatch("test/output/testGTF.gtf", gtfFile);

		// Test Interval13 model.
		gtfFile = "test/working/shortTrainInterval13.gtf";
		List<? extends TrainingSequence<?>> data3 = ih.readTrainingData("test/input/interval13/data/shortTrain.interval13.txt");
		oh = new OutputHandlerGeneCallStats(c.getModel(), ih);
		oh.writeGTF(data3, gtfFile);
		assertFilesMatch("test/output/shortTrain.gtf", gtfFile);	
	}
}


/** Below is an explanation of the fields in a GTF file:
 * 		1.  SEQNAME - The name of the sequence.  Typically a chomosome or a contig.
 * 		2.  SOURCE  - The program that generated this feature.
 * 		3.  FEATURE - The name of this type of feature.  I.e. "CDS", "start_codon", "stop_codon".
 * 		4.  START   - The starting position of the feature in the sequence.  The first base is 1.
 * 		5.  END     - The ending position of the feature (inclusive).
 * 		6.  SCORE   - A score between 0 and 1000.  If ther is no score value, enter ".".
 * 		7.  STRAND  - Valid entries include "+", "-", "." (for don't know).
 * 		8.  FRAME   - A number between 0-2 (inclusive) that represents the reading frame of the
 * 					  first base.
 * 		9.  GROUPING ATTRIBUTES - Attribute keys and values.
 * 
 *  More information on frames:
 *    Frame is the number of bases in this region before you get in frame.  That is, if it is 0,
 *    the first three bases in this element are a codon.  If it's 1, the first base is the end 
 *    of a codon hanging over from the last exon, and the next three are the first codon.  If 
 *    it's 2, the first two bases are the end of the previous codon, and the next three are the 
 *    first codon in this feature.  The first exon in each + stranded transcript has a frame of 0
 *    and the rest can vary all over the place.
 *    
 *  Note:  In this test the hidden sequence is used to generate the GTF file, the nucleotide sequence
 *    is not used.
 */

