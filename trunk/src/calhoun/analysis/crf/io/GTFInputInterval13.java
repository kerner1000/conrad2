package calhoun.analysis.crf.io;

import java.io.IOException;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang.builder.CompareToBuilder;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import calhoun.analysis.crf.FeatureList;
import calhoun.analysis.crf.features.interval13.GeneConstraintsInterval13;
import calhoun.analysis.crf.solver.CacheProcessor.FeatureEvaluation;
import calhoun.util.Assert;
import calhoun.util.CheckException;
import calhoun.util.FileUtil;

public class GTFInputInterval13 implements TrainingSequenceIO {
	private static final long serialVersionUID = 4413724139445660883L;
	private static final Logger LOGGER = LoggerFactory.getLogger(GTFInputInterval13.class);
	private static final Log log = LogFactory.getLog(GTFInputInterval13.class);
	
	String nameComponent = "name";
	
	static class CDS implements Comparable<CDS> {
		String seq;
		String gene;
		long start;
		long stop;
		char strand;
		
		public int compareTo(CDS other) {
			// Compare order is seq, gene, start.  Using gene second ensures we don't have duplicated gene names or overlapping genes
		     return new CompareToBuilder().append(seq, other.seq).append(gene, other.gene).append(start, other.start).toComparison();
		}

		public boolean equals(CDS other) {
		     return new EqualsBuilder().append(seq, other.seq).append(gene, other.gene).append(start, other.start).isEquals();
		}

		public long hashcode() {
		     return new HashCodeBuilder().append(seq).append(gene).append(start).hashCode();
		}
	}
	
	public void readTrainingSequences(Object location, List<TrainingSequence<Map<String, Object>>> seqs) throws IOException {
		String[][] gtf = FileUtil.readFlatFile((String) location);
	
		// Create a list of exons for each sequence
		Map<String, List<CDS>> exonLists = new HashMap();
		for(TrainingSequence<Map<String, Object>> seq : seqs) {
			exonLists.put((String) seq.getInputSequence().getComponent("name").getX(0), new ArrayList());
		}
		long threwout = 0;
		long kept = 0;
		// First read in all of the exons into an array and sort by sequence and position
		for(String[] row : gtf) {
			log.debug("current row: " + Arrays.asList(row));
			if(row[2].equalsIgnoreCase("cds")) {
				CDS cds = new CDS();
				cds.seq = row[0];
				cds.start = Integer.parseInt(row[3]);
				cds.stop = Integer.parseInt(row[4]);
				cds.strand = row[6].charAt(0);
				
				// Parse out the gene identifier
				String attributes = row[8];
				//log.warn(attributes);
				for(String pair : attributes.split(";")) {
					String[] keyValue = pair.trim().split("[ =]");
					
					// Strip quotes if they surround the ids
					try{
					if(keyValue[1].charAt(0)=='"')
						keyValue[1] = keyValue[1].substring(1,keyValue[1].length()-1);
					}catch(Exception e){
						e.printStackTrace();
						System.out.println(attributes);
					}
					// Check for something that indicates where this CDS belongs
					if(keyValue[0].equals("gene_id") || keyValue[0].equals("Parent")) {
						cds.gene = keyValue[1];
					}
					//log.warn("Key="+keyValue[0]+" Value="+keyValue[1]);
				}
				Assert.a(cds.gene != null);
				if (exonLists.containsKey(cds.seq)) {
					exonLists.get(cds.seq).add(cds);
					kept++;
				} else {
					threwout++;
				}
				
			}
		}
		if (threwout > 0) {
			log.warn("Threw out " + threwout + " of " + (threwout + kept) + " exons");
		}
		// Now go through and populate int vectors for each sequence 
		for(TrainingSequence<Map<String, Object>> seq : seqs) {
			String name = (String) seq.getInputSequence().getComponent("name").getX(0);
			List<CDS> exons = exonLists.get(name);
			
			// Sort in position order
			Collections.sort(exons);
			
			int[] states = new int[seq.length()];
			//log.warn("Seq: "+name+" Length: "+seq.length());
		
			
			mapExonsToStates(exons, states);
			seq.setY(states);

			//confirmSeq(seq);
		}
	}

	void mapExonsToStates(List<CDS> exons, int[] states) {
		if(exons.size() == 0)
			return;

		// 1-based index of the last base of the previous exon (or 0-based index of the first base of gap)
		long currentPosition = Long.MIN_VALUE;
		long exonState = Long.MIN_VALUE;
		long intronState = Long.MIN_VALUE;
		
		String currentGene = null;
		for(CDS exon : exons) {
			try{
//			log.warn(exon.gene+": "+exon.start+"-"+exon.stop + " currentPos="+currentPosition);
			Assert.a(exon.start > currentPosition);
			boolean sameGene = exon.gene.equals(currentGene); 
			if(sameGene) {
				// Gap was intron, fill in the previous state
				// Intergenic is the default 0, and so we don't fill that in.
				if(exon.strand == '+') {
					intronState = (3-(currentPosition - (exonState-1))%3)%3+4;
				}
				else {
					intronState = (currentPosition - (exonState-7))%3+10;
				}
				log.warn(String.format("%d-%d State: %d", currentPosition, exon.start -1, intronState));
				if(currentPosition < Integer.MIN_VALUE || currentPosition > Integer.MAX_VALUE)
					log.error("error: currentPosition out of range!");
				if(exon.start -1 < Integer.MIN_VALUE || exon.start -1 > Integer.MAX_VALUE)
					log.error("error: exon.start-1 out of range!");
				if(intronState < Integer.MIN_VALUE || intronState > Integer.MAX_VALUE)
					log.error("error: intronState out of range!");
				Arrays.fill(states, (int)currentPosition, (int)exon.start -1, (int)intronState);
			}
			// Fill in the current exon
			if(sameGene) {
				if(exon.strand == '+')
					exonState = (exon.start-1+intronState-4)%3 + 1;
				else
					exonState = ((exon.start-1) - (intronState-10))%3+7;
			}
			else {
				// New gene, only the current start matters
				exonState = ((exon.start-1)%3) + 1 + (exon.strand == '-' ? 6:0);
			}
			//log.warn(String.format("%d-%d State: %d", exon.start-1, exon.stop-1, exonState));
			if(exon.start-1 < Integer.MIN_VALUE || exon.start-1 > Integer.MAX_VALUE)
				log.error("error: exon.start-1 out of range!");
			if(exon.stop < Integer.MIN_VALUE ||exon.stop > Integer.MAX_VALUE)
				log.error("error: exon.stop out of range!");
			if(exonState < Integer.MIN_VALUE || exonState > Integer.MAX_VALUE)
				log.error("error: intronState out of range!");
			Arrays.fill(states, (int)exon.start-1, (int)exon.stop, (int)exonState);
			
			currentGene = exon.gene;
			currentPosition = exon.stop;
			//debug("currentPosition="+currentPosition);
			
			}catch(CheckException e){
				// ignore gene
			}
		}
	}
	
	public void writeTrainingSequences(Object location, Iterator<int[]> data) throws IOException {
	}

	/**
	 * @return Returns the nameComponent.
	 */
	public String getNameComponent() {
		return nameComponent;
	}

	/**
	 * @param nameComponent The nameComponent to set.
	 */
	public void setNameComponent(String nameComponent) {
		this.nameComponent = nameComponent;
	}

	/* This is debugging code that lets you get a better idea if problems occur.  Specific to interval13 */
	void confirmSeq(TrainingSequence<?> seq) {
		int[] states = seq.y;
		DirectFeatureList f = new DirectFeatureList();
		GeneConstraintsInterval13 g = new GeneConstraintsInterval13();
		InputSequenceCharacter a = (InputSequenceCharacter) seq.getInputSequence().getComponent("ref");
		for(int i = 1; i<seq.length(); ++i) {
			g.evaluateEdge(a, i, states[i-1], states[i], f);
			Assert.a(f.valid, String.format("Invalid at %d: %d-%d",i, states[i-1], states[i]));
		}
	}
	class DirectFeatureList implements FeatureList {
		FeatureEvaluation evals1;
		public int position;
		boolean valid = true;
		
		public DirectFeatureList() {
		}
		
		public void addFeature(int index, double val) {
			evals1.index[position] = (short) index;
			evals1.value[position++] = (float) val;
		}

		/** Returns the invalid flag. */
		public boolean isValid() {
			return valid;
		}

		/** Invalidates results. */
		public void invalidate() {
			valid = false;
		}
	}
	
	private static void info(Object msg){
		if(LOGGER.isInfoEnabled())
			LOGGER.info(msg.toString());
	}
	
	private static void debug(Object msg){
		if(LOGGER.isDebugEnabled())
			LOGGER.debug(msg.toString());
	}
	
	private static void warn(Object msg){
		LOGGER.warn(msg.toString());
	}
}
