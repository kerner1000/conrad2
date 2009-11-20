package calhoun.analysis.crf.features.supporting.phylogenetic;
import java.io.IOException;
import java.io.Serializable;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.util.Assert;

public class PhylogeneticTreeFelsensteinOrder implements Serializable {
	private static final long serialVersionUID = -2394371116379001213L;
	private static final Log log = LogFactory.getLog(PhylogeneticTreeFelsensteinOrder.class);
/* For now, hardcoded for the phylogeny and relative branch lengths of the
 * aligned cryptococcus serotypes.  However, all that is needed to make this
 * general is a parser for NSF formatted descriptions of trees and a comparison
 * to the order of species in the multiple alignment.
 */
	
	int[]     ileft,iright;
	double[]  bleft,bright;
	int nSpecies,nSteps,nNodes;
	
	public PhylogeneticTreeFelsensteinOrder() {
		// Below is a representation of the crypto tree which I got from Dan:
		// ((cnBB: 0.034084, cnBV: 0.034065): 0.186256, cnA: 0.081165, cnDT: 0.071992):0.0;
		// Below is the one I'll be using since it has all 5 species:
		// (((cnDT:0.0025,cnDS:0.0025):0.0475,cnAB:0.05):0.02,(cnBB:0.02,cnBV:0.02):0.05):0;

		// the order of species in the multiple alignment is cnDT,cnDS,cnA,cnBB,cnBV

		// Would be nice to figure out the sequence of four steps below given the string
		// representing the phylogenetic tree (NSF format) and the order of species in the multiple alignment column
		// 5 = (0,0.0025 ; 1,0.0025)
		// 6 = (5,0.0475 ; 2,0.05)
		// 7 = (3,0.02 ; 4,0.02)
		// 8 = (6,0.02 ; 7,0.05)

		ileft = new int[]{0,5,3,6};
		iright = new int[]{1,2,4,7};
		bleft = new double[]{0.0025,0.0475,0.02,0.02};
		bright = new double[]{0.0025,0.05,0.02,0.05};
		
		setup();
	}

	public PhylogeneticTreeFelsensteinOrder(int[] ileft, int[] iright, double[] bleft, double[] bright) {
		this.ileft = ileft;
		this.iright = iright;
		this.bleft = bleft;
		this.bright = bright;
		
		setup();
	}
	
	public void summarize() throws IOException {
		log.debug("Summary of an order-of-computations object for Felsenstein recursion");
		log.debug("Number of species: " + nSpecies);
		log.debug("Number of steps:   " + nSteps);
		log.debug("Number of nodes:   " + nNodes);
		log.debug("  Nodes 0 to " + (nSpecies-1) + " are given.");
		//for (int j=0; j<nSpecies; j++) {
		//	stream.println("  Node " + j + " is given");
		//}
		for (int step=0; step<nSteps; step++) {
			log.debug("  Node " + (nSpecies + step) + " combines nodes " + ileft[step] + " and " + iright[step] +
					" with branchlengths " + bleft[step] + " and " + bright[step]);
		}
	}
		
	private void setup() {
		nSteps = ileft.length;
		Assert.a(iright.length == nSteps);
		Assert.a(bleft.length == nSteps);
		Assert.a(bright.length == nSteps);
		
		nSpecies = nSteps+1;
		nNodes = nSpecies + nSteps;
	}
	
	public int numSpecies() {
		return nSpecies;
	}

	public int numSteps() {
		return nSteps;
	}
	
	public int numNodes() {
		return nNodes;
	}
	
	public int[] getileft() {
		return ileft;
	}
	
	public int[] getiright() {
		return iright;
	}

	public double[] getbleft() {
		return bleft;
	}

	public double[] getbright() {
		return bright;
	}
	
}

