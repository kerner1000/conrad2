package calhoun.analysis.crf.test;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.features.supporting.phylogenetic.PhylogeneticTreeFelsensteinOrder;
import calhoun.analysis.crf.features.supporting.phylogenetic.RootedBinaryPhylogeneticTree;
import calhoun.util.AbstractTestCase;

public class NewickTreeManipulationTest extends AbstractTestCase {
	private static final Log log = LogFactory.getLog(CRFIOTest.class);
	boolean debug = log.isDebugEnabled();
	
	public void testAspTree() throws Exception {
		String treeStr = "(((AC2:1.0,(NF2:1.0,AF2:1.0):1.0):1.0,((AFL2:1.0,A_oryzae_RIB40:1.0):1.0,AT1:1.0):1.0):1.0,AN1:1.0):0.0;";
		RootedBinaryPhylogeneticTree tree = new RootedBinaryPhylogeneticTree(treeStr);	
		assertEquals(7, tree.getNumSpecies());
	}
	
	public void testBigTree() throws Exception {
		RootedBinaryPhylogeneticTree RBPT = ExampleTrees.bigTree();	
		RBPT.summarize_tree();
		
		String[]  selected = new String[]{"elephant","human","mouse","dog","snufalupagus"};
		
		RootedBinaryPhylogeneticTree subtree = RBPT.subtree(selected);
		subtree.summarize_tree();	
	}
	
	public void testCryptoTree() throws Exception {
		RootedBinaryPhylogeneticTree RBPT = ExampleTrees.crypto5();	
		RBPT.summarize_tree();
		
		String[] msaOrder = new String[]{"cnDT","cnDS","cnAB","cnBB","cnBV"};
		PhylogeneticTreeFelsensteinOrder fo = RBPT.getFelsensteinOrder(msaOrder);
		
		fo.summarize();
		
	}
	
	public void testEncodeSpeciesFromBigTree() throws Exception {
		RootedBinaryPhylogeneticTree RBPT = ExampleTrees.bigTreeForEncode();	
		RBPT.summarize_tree();

		String[] speciesList = new String[]{"hg17","hedgehog","shrew","baboon","tenrec",
				"owl_monkey","rfbat","dusky_titi","chimp","zebrafish","platypus","cow",
				"mouse","macaque","chicken","elephant","colobus_monkey","marmoset",
				"mouse_lemur","armadillo","rat","tetraodon","rabbit","monodelphis",
				"galago","xenopus","dog","fugu"};
		
		RootedBinaryPhylogeneticTree subtree = RBPT.subtree(speciesList);
		subtree.summarize_tree();	
		
		System.out.println("Now computing a Felsenstein ordering for the subtree:");
		subtree.getFelsensteinOrder(speciesList).summarize();
		
	}
	
}
