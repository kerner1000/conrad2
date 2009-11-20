package calhoun.analysis.crf.io;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.io.IntervalInputSequence.IntervalRangeMapValue;
import calhoun.analysis.crf.statistics.PredictedActualBinaryContingencyTable;
import calhoun.util.Assert;
import calhoun.util.DenseBooleanMatrix2D;
import calhoun.util.FileUtil;
import calhoun.util.RangeMap;

/** A legacy output handler that computes basic stats, gene calling statsm and then writes out a GTF file.
 */
public class OutputHandlerGeneCallPredict implements OutputHandler {
	private static final long serialVersionUID = 2014487490985409134L;

	private static final Log log = LogFactory.getLog(OutputHandlerGeneCallPredict.class);

	private ModelManager manager;
	private InputHandler inputHandler;
	String location;
	boolean writeTrainingData = false;

	/** default constructor.  <code>ModelManager</code> and <code>InputHandler</code> must be configured separately. */
	public OutputHandlerGeneCallPredict() {
	}
	
	/** creates an output handler using this model and input handler
	 * @param manager the model used for gene calling.  Used for calculating stats.
	 * @param inputHandler input handler which will be used for writing out the input sequence with the results.
	 */
	public OutputHandlerGeneCallPredict(ModelManager manager, InputHandler inputHandler) {
		this.inputHandler = inputHandler;
		setManager(manager);
	}
	
	public void setOutputLocation(String location) {
		this.location = location;
	}
	
	public void writeOutput(InputSequence<?> sequence, int[] hiddenStates) throws IOException {
		throw new UnsupportedOperationException();
	}

	public void writeTestOutput(InputSequence<?> sequence, int[] truePath, int[] hiddenStates) throws IOException {
		calcResultIncrement(new TrainingSequence(sequence, truePath), hiddenStates);
	}

	public void outputComplete() throws IOException {
		if(location != null) {
			if(writeTrainingData) {
				try {
					inputHandler.writeTrainingData(location, labeled);
				}
				catch(Exception ex) {
					log.warn("Unable to write training data", ex);
				}
			}
			writeGTF(labeled, location + ".gtf");
			System.out.print(this);
			writeResults(location + ".dat");
		}
	}

	/** retursn the input handler used to write out the input sequences
	 * @return the inputHandler which will be used to write out the input sequences
	 */
	public InputHandler getInputHandler() {
		return inputHandler;
	}

	/** sets the inputHandler used to write out the input sequences
	 * @param inputHandler the inputHandler used to write out the input sequences
	 */
	public void setInputHandler(InputHandler inputHandler) {
		this.inputHandler = inputHandler;
	}

	/** gets the model used to generate results
	 * @return the model used to generate results
	 */
	public ModelManager getManager() {
		return manager;
	}

	/** sets the model used to generate results
	 * @param manager the model used to generate results
	 */
	public void setManager(ModelManager manager) {
		this.manager = manager;

		ctCodingNucleotide = new PredictedActualBinaryContingencyTable();
			
		ctExons = new PredictedActualBinaryContingencyTable();
		ctExons.forgetTN();
		
		nStates = manager.getNumStates();
		ctStates = new ArrayList<PredictedActualBinaryContingencyTable>();
		for (int i=0; i<nStates; i++) {
			ctStates.add(new PredictedActualBinaryContingencyTable());
		}

		DenseBooleanMatrix2D LT = manager.getLegalTransitions();
		fromInd = new ArrayList<Integer>();
		toInd   = new ArrayList<Integer>();
		for (int from=0; from<nStates; from++) {
			for (int to=0; to<nStates; to++ ) {
				if (LT.getQuick(from,to)) {
					fromInd.add(from);
					toInd.add(to);
				}
			}
		}
		nTransitions = fromInd.size();
		ctTransitions = new ArrayList<PredictedActualBinaryContingencyTable>();
		for (int i=0; i<nTransitions; i++) {
			ctTransitions.add(new PredictedActualBinaryContingencyTable());
		}
	}
			
	private List<TrainingSequence<?>> labeled = new ArrayList<TrainingSequence<?>>();
	private int correct = 0;
	private int incorrect = 0;
	private int perfect = 0;
	private int imperfect = 0;
	private transient double[] viterbiScores;
	private double lla = 0, llv = 0;

	// Info we'll need to know about the model in order to be clever about gathering stats
	// (the manager must be provided, other things derived from it)
	private int nStates;
	private int nTransitions;
	private List<Integer> fromInd;
	private List<Integer> toInd;;
	
	// The 2x2 contingency tables for which we'll keep track of results:
	private PredictedActualBinaryContingencyTable ctCodingNucleotide;
	private PredictedActualBinaryContingencyTable ctExons;	
	private List<PredictedActualBinaryContingencyTable> ctStates;
	private List<PredictedActualBinaryContingencyTable> ctTransitions;
	
	/** returns the exact nucleotide accuracy of the result */
	public float getAccuracy() {
		return correct / (float)(correct+incorrect);
	}

	public static class Results implements Serializable { 
		private static final long serialVersionUID = 9082449588200635355L;
		public PredictedActualBinaryContingencyTable ctCodingNucleotide;
		public PredictedActualBinaryContingencyTable ctExons;	
		public List<PredictedActualBinaryContingencyTable> ctStates;
		public List<PredictedActualBinaryContingencyTable> ctTransitions;
		public int correct;
		public int incorrect;
		public int perfect;
		public int imperfect;
	}
	
	void writeResults(String loc) throws IOException {
		Results results = new Results();
		results.ctCodingNucleotide = ctCodingNucleotide;
		results.ctExons = ctExons;
		results.ctStates = ctStates;
		results.ctTransitions = ctTransitions;
		results.correct = correct;
		results.incorrect = incorrect;
		results.perfect = perfect;
		results.imperfect = imperfect;
		FileUtil.writeObject(loc, results);
	}
	
	@Override
	public String toString() {
		String ret = "";
		
		for (int s=0; s<nStates; s++) {
			ret += "[State=" + manager.getStateName(s) + "] ";
			ctStates.get(s).freeze();
			ret += "Predicted: " + ctStates.get(s).pp();
			ret += "\n";
		}
		
		for (int t=0; t<nTransitions; t++) {
			ret += "[Transition " + manager.getStateName(fromInd.get(t)) + " --> " + manager.getStateName(toInd.get(t)) + " ] ";
			ctTransitions.get(t).freeze();
			ret += "Predicted: " + ctTransitions.get(t).pp();
			ret += "\n";
		}		

		ctCodingNucleotide.freeze();
		ret += "[Coding nucleotides] Predicted: " + ctCodingNucleotide.pp() + "\n";
		
		ctExons.freeze();
		ret += "[Coding exons] Predicted: " + ctExons.pp() + "\n";		
		
		if (lla>0) {
			ret += "LLA:" + lla + "  LLV:" + llv + "  " + "\n";
		}
		
		//ret += String.format("Perfectly predicted hidden sequences: %d/%d %.2f %%",perfect,perfect+imperfect,perfect*100.0/(float) (perfect+imperfect))+ "\n";
		
		//ret += String.format("Nucleotide Hidden State Agreement: %d/%d %.2f %%",correct, correct + incorrect, correct * 100.0 / (float) (correct + incorrect)) + "\n";

		
		return ret;
	}
	
	/** calculates statstics and output for results on a given test sequence */
	public void calcResultIncrement(TrainingSequence training, int[] predictedHiddenSequence) {
		labeled.add(new TrainingSequence(training.getInputSequence(), predictedHiddenSequence));  // This is only place that labelled gets added to???
		// So I guess the results just get built up incrementally, both the actuall hidden sequences and the stats?
		Assert.a(training.length() == predictedHiddenSequence.length);
		int[] actualHiddenSequence = new int[training.length()];
		for (int i=0; i<training.length(); i++) {
			actualHiddenSequence[i] = training.getY(i);
		}
		boolean thisperfect = true;
		for (int i = 0; i < predictedHiddenSequence.length; ++i) {
			int predY = predictedHiddenSequence[i];
			int realY = actualHiddenSequence[i];

			if (realY == predY) {	 correct += 1; } else { incorrect += 1; thisperfect = false; }
			
			ctCodingNucleotide.increment(isCodingPlus(predY),isCodingPlus(realY));
			ctCodingNucleotide.increment(isCodingMinus(predY),isCodingMinus(realY));
			
			for (int s=0; s<nStates; s++) {
				ctStates.get(s).increment((predY==s),(realY==s));
			}	
		}
		if (thisperfect) {
			perfect++;
		} else {
			imperfect++;
		}
		for (int i = 1; i < predictedHiddenSequence.length; ++i) {
			int predY = predictedHiddenSequence[i];
			int realY = actualHiddenSequence[i];
			int predYp = predictedHiddenSequence[i-1];
			int realYp = actualHiddenSequence[i-1];			

			for (int t=0; t<nTransitions; t++) {
				boolean bPred = ( (predYp==fromInd.get(t)) && (predY==toInd.get(t)) );
				boolean bReal = ( (realYp==fromInd.get(t)) && (realY==toInd.get(t)) );
				ctTransitions.get(t).increment(  bPred  ,  bReal  );
			}	
		}

		// Now let's increment the contingency table for exons; note that here not counting TN's
		RangeMap predExonsPlus = new RangeMap();
		RangeMap predExonsMinus = new RangeMap();
		RangeMap realExonsPlus = new RangeMap();
		RangeMap realExonsMinus = new RangeMap();
		makeExonRangeMapFrom13SV(predictedHiddenSequence,predExonsPlus,predExonsMinus);
		makeExonRangeMapFrom13SV(actualHiddenSequence,realExonsPlus,realExonsMinus);
		incrementCTFromRangeMaps(ctExons,predExonsPlus,realExonsPlus);
		incrementCTFromRangeMaps(ctExons,predExonsMinus,realExonsMinus);
	
	}

	private void incrementCTFromRangeMaps(PredictedActualBinaryContingencyTable ct, RangeMap pred, RangeMap real) {
		// By looping through the predictions, can get at TP and FP
		Set<IntervalRangeMapValue> pv = pred.values();
		Iterator<IntervalRangeMapValue> pvi = pv.iterator();
		while(pvi.hasNext()) {
			IntervalRangeMapValue irmv = pvi.next();
			Set<IntervalRangeMapValue> vals = real.find(irmv.start,irmv.end);
			if(vals.size() == 0) {
				ct.incrementFP();
			}
			else {
				IntervalRangeMapValue val = vals.iterator().next();
				if(val.start == irmv.start && val.end == irmv.end) {
					ct.incrementTP();
				} else {
					ct.incrementFP();
				}
			}
		}
		Set<IntervalRangeMapValue> rv = real.values();
		Iterator<IntervalRangeMapValue> rvi = rv.iterator();
		while(rvi.hasNext()) {
			IntervalRangeMapValue irmv = rvi.next();
			if (!pred.hasEntry(irmv.start,irmv.end)) {
				ct.incrementFN();
			}
			Set<IntervalRangeMapValue> vals = pred.find(irmv.start,irmv.end);
			if(vals.size() == 0) {
				ct.incrementFN();
			}
			else {
				IntervalRangeMapValue val = vals.iterator().next();
				if(val.start == irmv.start && val.end == irmv.end) {
				} else {
					ct.incrementFN();
				}
			}
		}	
	}


	private void makeExonRangeMapFrom13SV(int[] hidden, RangeMap exonsPlus, RangeMap exonsMinus) {
		
		int len = hidden.length;
		
		for (int i=1; i<len; i++) {
			if ((!isCodingPlus(hidden[i-1]) && (isCodingPlus(hidden[i])))) {
				int j=i;
				while ((isCodingPlus(hidden[j])) &&(j<(len-1))) { j++; }
				exonsPlus.add(i,j,new IntervalRangeMapValue(i,j,1.0));
				//log.info("Add + "+i+" "+j);
			}
			if ((!isCodingMinus(hidden[i-1]) && (isCodingMinus(hidden[i])))) {
				int j=i;
				while ((isCodingMinus(hidden[j])) &&(j<len-1)) { j++; }
				exonsMinus.add(i,j,new IntervalRangeMapValue(i,j,1.0));
				//log.info("Add - "+i+" "+j);
			}
		}
	}

	private boolean isCodingPlus(int y) {
		Assert.a( (y>=0) && (y<13) );
		if ( (y==1) || (y==2) || (y==3) ) { return true; }
		return false;
	}

	private boolean isCodingMinus(int y) {
		Assert.a( (y>=0) && (y<13) );
		if ( (y==7) || (y==8) || (y==9) ) { return true; }
		return false;
	}
	
	public void loglikelihoodIncrement(double logLikelihoodActual, double logLikelihoodViterbi) {
		lla += logLikelihoodActual;
		llv += logLikelihoodViterbi;
	}

	public TrainingSequence getLabeled(int i) {
		return labeled.get(i);
	}
	
	String seqName;
	String genePrefix;
	long   offset;
	
	// This function converts a 13 state model hidden sequence to a GTF file.  
	public void writeGTF(List<? extends TrainingSequence<?>> refStates, String filename) throws IOException
	{		
		int ref, geneNum, seqCount, frame=-1;
		long i, exonStart, exonEnd, end;
		boolean inPlusExon, inMinusExon, firstExon, startCodonSplit;
		String strand;
		Writer fout = new BufferedWriter(new FileWriter(filename));	
		exonStart = exonEnd = 0;
		geneNum = 1;
		seqCount = 0;
		
		// Determine if model is tricycle13 or interval13.
		boolean interval13 = false;
		int prevState, state;
		for (TrainingSequence seq : refStates) {
			if (seq.length() == 0) 	 continue;
			
			prevState = seq.getY(0);
			for (i=1; i<seq.length(); i++) {
				state = seq.getY((int)i);
				if (prevState == 0 && (state==2 || state==3 || state==7 || state==8)) {
					interval13 = true;
					break;
				}
				prevState = state;
			}
			if (interval13)
				break;
		}
		
		for (TrainingSequence seq : refStates) {	
			
			if (interval13) {
				SequenceConverter.convertSeqFromInterval13ToTricycle13(seq);
			}
			
			inPlusExon  = false;
			inMinusExon = false;
			firstExon   = true;
			startCodonSplit = false;
			
			parseSeqName(seq, seqCount);
									
			for (i=0; i<seq.length(); i++)
			{			
				ref = seq.getY((int)i);
				
				if (ref == 1 || ref == 2 || ref == 3)		// in a plus exon
				{
					if (!inPlusExon)
					{
						exonStart = i+1;
						inPlusExon = true;
						frame = setFrame(ref);
					}
				}
				else if (ref == 7 || ref == 8 || ref == 9)	// in a minus exon
				{
					if (!inMinusExon)
					{
						exonStart = i+1;
						inMinusExon = true;
						frame = setFrame(ref);
						if (firstExon) {
							if (i < 3)
								System.err.println("Minus strand gene start is within 3 nucleotides of sequence start.  No stop codon writen to GTF for gene starting at position " + (exonStart+offset));
							else
								writeGFTLine(fout,seqName,"stop_codon",exonStart+offset-3,exonStart+offset-1,"-",frame,genePrefix,geneNum);						
	
						}
					}
				}
				else if ( inPlusExon  && (ref == 4  || ref == 5  || ref == 6) )	{ // just ended an exon on plus strand, now in a plus intron
					strand = "+";	
					inPlusExon = false;
					exonEnd = i;
					if (firstExon) {
						if (exonEnd - exonStart + 1 < 3)	{ end = exonEnd + offset; startCodonSplit = true;}
						else								{ end = exonStart+offset+2; }
						writeGFTLine(fout,seqName,"start_codon",exonStart+offset,end,strand,frame,genePrefix,geneNum);						
						firstExon = false;
					}
					else if (startCodonSplit) {	// at second exon that contains part of start codon
						Assert.a(frame==1 || frame==2);
						writeGFTLine(fout,seqName,"start_codon",exonStart+offset,exonStart+offset+frame-1,strand,frame,genePrefix,geneNum);						
						startCodonSplit = false;
					}
					writeGFTLine(fout,seqName,"CDS",exonStart+offset,exonEnd+offset,strand,frame,genePrefix,geneNum);
				}
				else if (inMinusExon && (ref == 10 || ref == 11 || ref == 12))  { // just ended an exon on minus strand, now in a minus intron
					strand = "-";
					inMinusExon = false;
					firstExon = false;
					exonEnd = i;
					writeGFTLine(fout,seqName,"CDS",exonStart+offset,exonEnd+offset,strand,frame,genePrefix,geneNum);
				}
				else								// now in intergenic region
				{
					boolean write = true;
					if (inPlusExon)					// was in gene at previous nucleotide
					{
						strand = "+";
						exonEnd = i;			
						if (firstExon) {
							if (exonEnd - exonStart + 1 < 3) {
								System.err.println("Single '" + strand + "' strand exon is < 3 bases for sequence '" + seqName + "'.  exonStart=" + exonStart + "  exonEnd=" + exonEnd);
								write = false;
							}
							else {
								writeGFTLine(fout,seqName,"start_codon",exonStart+offset,exonStart+offset+2,strand,frame,genePrefix,geneNum);
							}
						}
						else if (startCodonSplit) {	// at second exon that contains part of start codon
							Assert.a(frame==1 || frame==2);
							writeGFTLine(fout,seqName,"start_codon",exonStart+offset,exonStart+offset+frame-1,strand,frame,genePrefix,geneNum);						
						}
						if (write) {
							writeGFTLine(fout,seqName,"CDS",exonStart+offset,exonEnd+offset,  strand,frame,genePrefix,geneNum);
							writeGFTLine(fout,seqName,"stop_codon",exonEnd+offset+1,exonEnd+offset+3,strand,0,    genePrefix,geneNum);
						}
						inPlusExon  = false;
						firstExon   = true;
						startCodonSplit = false;
						geneNum++;
					}
					else if (inMinusExon) {
						strand = "-";
						long prevExonEnd = exonEnd;
						exonEnd = i;			
						if (firstExon && exonEnd - exonStart + 1 < 3) {
							System.err.println("Single '" + strand + "' strand exon is < 3 bases for sequence '" + seqName + "'.  exonStart=" + exonStart + "  exonEnd=" + exonEnd);
						}
						else if (exonEnd - exonStart + 1 < 3) {	// this exon is < 3 bases, need to split start codon
							if (exonEnd - exonStart + 1 == 2) { // this exon is 2 bases
								writeGFTLine(fout,seqName,"start_codon",prevExonEnd+offset,prevExonEnd+offset,strand,0,genePrefix,geneNum);				
								writeGFTLine(fout,seqName,"CDS",exonStart+offset,exonEnd+offset,strand,frame,genePrefix,geneNum);
								writeGFTLine(fout,seqName,"start_codon",exonStart+offset,exonEnd+offset,strand,2,genePrefix,geneNum);				
							}
							else if (exonEnd - exonStart + 1 == 1) 	{ // this exon is 1 base
								writeGFTLine(fout,seqName,"start_codon",prevExonEnd+offset-1,prevExonEnd+offset,strand,0,genePrefix,geneNum);				
								writeGFTLine(fout,seqName,"CDS",exonStart+offset,exonEnd+offset,strand,frame,genePrefix,geneNum);
								writeGFTLine(fout,seqName,"start_codon",exonStart+offset,exonEnd+offset,strand,1,genePrefix,geneNum);												
							}
						}
						else {
							writeGFTLine(fout,seqName,"CDS",exonStart+offset,exonEnd+offset,strand,frame,genePrefix,geneNum);
							writeGFTLine(fout,seqName,"start_codon",exonEnd+offset-2,exonEnd+offset,strand,0,    genePrefix,geneNum);				
						}			
						inMinusExon = false;
						firstExon   = true;
						startCodonSplit = false;
						geneNum++;			
					}
				}
			}
			seqCount++;
		}
		fout.close();
	}
	
	private void parseSeqName(TrainingSequence seq, int seqNum) {
		NameInputSequence nameInput = null;

		InputSequence<?> inputSeq = seq.getInputSequence();
		if(inputSeq instanceof InputSequenceComposite) {
			nameInput = (NameInputSequence) inputSeq.getComponent("name");
			
		}
		if(nameInput == null) {
			log.debug("Sequence name not specified.  Setting sequence name to 'SEQ_" + String.valueOf(seqNum) + "'");
			seqName    = "SEQ_" + String.valueOf(seqNum);	// Create a name and return.
			genePrefix = "SEQ_" + String.valueOf(seqNum);
			offset = 0;
			return;
		}
		String name =  nameInput.getName().trim();
		
		int colon1, colon2, numColons;
		
		if (name.startsWith("group:") || name.startsWith("seq:") ) {
			numColons = numOccurrences(name, ':');
			if (numColons == 1) {
				colon1 = name.indexOf(":");
				seqName = name;
				genePrefix = name.substring(colon1 + 1, name.length());
				offset = 0;
				return;
			}
			else if (numColons == 2) {
				colon1 = name.indexOf(":");
				colon2 = name.lastIndexOf(":");
				seqName = name.substring(0, colon2);
				genePrefix = name.substring(colon1 + 1, colon2);
				int pound = genePrefix.indexOf("#");
				if (pound > 0) {
					genePrefix = genePrefix.substring(0, pound);
				}
				setOffset(name.substring(colon2+1, name.length()));
				return;
			}
		}
		log.debug("Sequence name is in unexpected format.  Setting offset=0 and sequence name='" + name + "'.");
		seqName    = name;
		genePrefix = name;
		offset = 0;
	}
	
	// Returns the number of times the character 'c' occurs in 'str'
	private static int numOccurrences(String str, char c) {
		int num = 0;
		int index = str.indexOf(c);
		while (index != -1) {
			num++;
			index = str.indexOf(c, index+1);
		}
		return num;
	}

	private void setOffset(String str) {
		int numDashes, dash;
		numDashes = numOccurrences(str, '-');
		
		if (numDashes == 0) {
			offset = 0;
		}
		else if (numDashes == 1) {
			try {
				dash = str.indexOf("-");
				offset = Long.valueOf(str.substring(0, dash)) - 1;
			}
			catch (NumberFormatException e) {
				System.err.println("Sequence range values in unexpected format.  Setting offset=0 for sequence='" + seqName + "'.");
				offset = 0;
			}
		}
		else {
			System.err.println("Sequence range values in unexpected format.  Setting offset=0 for sequence='" + seqName + "'.");
			offset = 0;
		}
	}
	
	// Frame is the nmber of bases in this region befor you get in frame.  
	// That is, if frame is 0, the first three bases in this element are a codon.
	// If frame is 1, the first base is the end of a codon hanging over from the 
	//     end of the previous codon and the next three are the first codon in this feature.
	// If frame is 2, the first two bases are the end of the previous codon and the 
	//     next three are the first codon in this feature.
	private static int setFrame(int ref) {
		int frame = -1;
		
		switch (ref) {
		case 1:  frame = 0;  break;
		case 2:  frame = 2;  break;
		case 3:  frame = 1;  break;
		case 7:  frame = 1;  break;
		case 8:  frame = 2;  break;
		case 9:  frame = 0;  break;
		default:  Assert.a(false, "Error setting frame, ref = ", ref);
		}
		return frame;
	}

	// Outputs one line to the GTF file.  
	// NOTE:  source is assumed to be 'CONRAD', and score is assumed to be unknown and set to '.'.
	private static void writeGFTLine(Writer out, String seqName, String feature, long exonStart, long exonEnd, 
									 String strand, int frame, String genePrefix, int geneNum) throws IOException {

		Assert.a(frame==0 || frame==1 || frame==2, "Frame value invalid, frame = ", frame);
		
		String geneId = genePrefix + "G_" + String.valueOf(geneNum);
		String transId = genePrefix + "T_" + String.valueOf(geneNum) + ".1";
		
		out.write(seqName + "\t" + "CONRAD" + "\t" + feature + "\t" + exonStart + "\t" + exonEnd + "\t" +
				  "." + "\t" + strand + "\t" + frame + "\t" + 
				  "gene_id \"" + geneId + "\"; transcript_id \"" + transId + "\";\n");	
	}

	public double[] getViterbiScores() {
		return viterbiScores;
	}

	/**
	 * @return Returns the writeTrainingData.
	 */
	public boolean isWriteTrainingData() {
		return writeTrainingData;
	}

	/**
	 * @param writeTrainingData The writeTrainingData to set.
	 */
	public void setWriteTrainingData(boolean writeTrainingData) {
		this.writeTrainingData = writeTrainingData;
	}

}
