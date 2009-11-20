package calhoun.analysis.crf.io;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import calhoun.analysis.crf.features.supporting.phylogenetic.RootedBinaryPhylogeneticTree;
import calhoun.util.Assert;
import calhoun.util.FileUtil;

public class AlignmentTree implements InputComponentIO {
	private static final long serialVersionUID = 760405914814389112L;

	String component;

	public List<String> getComponentNames() {
		List<String> ret = new ArrayList();
		ret.add(component);
		return ret;
	}

	public void readInputSequences(String location, List<Map<String, InputSequence<?>>> inputs) throws IOException {
		Assert.a(inputs.size() > 0, "AlignmentTree can't be the first input.");

		// Read the tree from the file
		String[] contents = FileUtil.readFile(location).split("\n");
		
		RootedBinaryPhylogeneticTree tree = new RootedBinaryPhylogeneticTree(contents[0]);

		// Add it into every input sequence
		for(Map<String, InputSequence<?>> input : inputs) {
			input.put(component, new MultipleAlignmentInputSequence(contents[1].trim(), tree));
		}
	}

	public void writeInputSequences(String location, List<? extends Map<String, ? extends InputSequence<?>>> inputComponents) throws IOException {
		MultipleAlignmentInputSequence alignment = (MultipleAlignmentInputSequence) inputComponents.get(0).get(component);
		String tree = alignment.getTree().newick();
		FileUtil.writeFile(location, tree);
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
