package calhoun.analysis.crf.io;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.features.supporting.phylogenetic.RootedBinaryPhylogeneticTree;
import calhoun.util.Assert;

/** an input component that reads in multiple alignment sequences.  Creates 
 * {@link MultipleAlignmentInputSequence} objects for each sequence.  Can be used as
 * a regular input or as part of an interleaved input file.
 */
public class MultipleAlignmentInput extends InterleavedInputComponentBase {
	private static final long serialVersionUID = 1796622784861590659L;
	private static final Log log = LogFactory.getLog(MultipleAlignmentInput.class);

	public boolean read(BufferedReader r, Map<String, InputSequence<?>> output) throws IOException {
		List<String> speciesNames = new ArrayList<String>();
		List<String> consensuses = new ArrayList<String>();
		int nSpecies;
		
		String temp = r.readLine();
		if(temp == null)
			return false;
		
		try {
			nSpecies = Integer.parseInt(temp);
		} catch (Exception e){
			log.error("Offending line was : " + temp);
			throw new IOException();
		}
		
		temp = r.readLine();
		RootedBinaryPhylogeneticTree tree = new RootedBinaryPhylogeneticTree(temp);
		
		log.debug("Number of species : " + nSpecies);	
		
		for (int spec=0; spec<nSpecies; spec++) { 
			String str1 = r.readLine();
			Assert.a(str1 != null);

			log.debug("One of the species: " + str1);
			speciesNames.add(str1.substring(1));

			String str2 = r.readLine();
			Assert.a(str2 != null);
			consensuses.add(str2);
		}

		InputSequence<?> inputSeq = new MultipleAlignmentInputSequence(speciesNames, consensuses, speciesNames.get(0), tree);
		output.put(name, inputSeq);
		return true;
	}
	

	public void write(Writer w, Map<String, ? extends InputSequence<?>> data) throws IOException {
		MultipleAlignmentInputSequence inputSeq = (MultipleAlignmentInputSequence) data.get(name);
		w.write("" + inputSeq.nSpecies + "\n");
		w.write("" + inputSeq.tree.newick() + "\n");
		for (int spec=0; spec<inputSeq.nSpecies; spec++) {
			w.write(">" + inputSeq.speciesNames.get(spec));
			w.write('\n');
			w.write(inputSeq.consensuses.get(spec));
			w.write('\n');			
		}
	}
}
