package calhoun.analysis.crf.io;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import calhoun.seq.FastaIterator;
import calhoun.seq.FastaSequence;
import calhoun.seq.FastaWriter;
import calhoun.util.Assert;

public class FastaAlignmentInput implements InputComponentIO {
	private static final long serialVersionUID = 760405914814389112L;

	String component;

	public List<String> getComponentNames() {
		return Collections.singletonList(component);
	}

	public void readInputSequences(String location, List<Map<String, InputSequence<?>>> inputs) throws IOException {
		FastaIterator it = new FastaIterator(location);
		int seqNum = 0;
		String[] current = parseSeq(it.next());
		
		List<String> species = null;
		
		for(Map<String, InputSequence<?>> input : inputs) {
			MultipleAlignmentInputSequence alignment = (MultipleAlignmentInputSequence) input.get(component);
			Assert.a(alignment != null, "AlignmentTree must be input before the alignment Fasta");

			if(species == null) {
				species = new ArrayList();
				species.addAll(alignment.getTree().getSpeciesSet());
				
			}
			List<String> consensuses = new ArrayList<String>();
			for(int i=0; i<alignment.getNumSpecies(); ++i) {
				consensuses.add(null);
			}
			
			// Loads all species for this sequence
			int len = -1;
			String thisSeq = current[0];
			while(thisSeq.equals(current[0])) {
				int ix = species.indexOf(current[1]);
				Assert.a(ix != -1, "Seq: "+current[0]+". Species is missing: "+current[1]);
				consensuses.set(ix, current[2]);
				len = current[2].length();
				if(it.hasNext())
					current = parseSeq(it.next());
				else
					break;
			}
			
			// Pad missing species with gaps
			Assert.a(len != -1, "No alignments available for this input.");
			char[] gapChars = new char[len];
			Arrays.fill(gapChars, '-');
			String gaps = new String(gapChars);
			for(int i = 0; i < consensuses.size(); ++i) {
				if(consensuses.get(i) == null)
					consensuses.set(i, gaps);
			}

			alignment.setSpeciesAndConsensuses(species, consensuses);
			seqNum += 1;
		}
	}

	String[] parseSeq(FastaSequence seq) {
		String[] ret = new String[3];
		String[] header = seq.getHeader().split(" ");
		ret[0] = header[0];
		ret[1] = header[1];
		ret[2] = seq.getSequence();
		return ret;
	}
	
	public void writeInputSequences(String location, List<? extends Map<String, ? extends InputSequence<?>>> inputComponents) throws IOException {
		FastaWriter w = new FastaWriter(location, false);
		for(Map<String, ? extends InputSequence<?>> input : inputComponents) {
			MultipleAlignmentInputSequence alignment = (MultipleAlignmentInputSequence) input.get(component);
			List<String> names = alignment.getSpeciesNames();
			List<String> consensuses = alignment.getConsensusSeqs();
			for(int i=0; i<names.size(); ++i) {
				w.writeSeq(names.get(i), consensuses.get(i));
			}
		}
		
		w.close();
	}

	/**
	 * @return Returns the header.
	 */
	public String getComponent() {
		return component;
	}

	/**
	 * @param header The header to set.
	 */
	public void setComponent(String header) {
		this.component = header;
	}
}
