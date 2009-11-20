package calhoun.analysis.crf.io;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.apache.commons.lang.StringUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.features.supporting.phylogenetic.PhylogeneticTreeFelsensteinOrder;
import calhoun.analysis.crf.features.supporting.phylogenetic.RootedBinaryPhylogeneticTree;
import calhoun.seq.KmerHasher;
import calhoun.util.Assert;

/** an input sequence where each element represents one column of a multiple alignment.
 */
public class MultipleAlignmentInputSequence implements InputSequence<MultipleAlignmentInputSequence.MultipleAlignmentColumn> {
	private static final Log log = LogFactory.getLog(MultipleAlignmentInputSequence.class);

	// Raw original data
	List<String> speciesNames;
	List<String> consensuses;
	RootedBinaryPhylogeneticTree tree;
	String refSpecies;
	
	// Derived data
	int nSpecies;
	int consensusLength;
	KmerHasher h;
	ArrayList<Integer> ref2con;
	int[] con2refLeft;
	int[] con2refRight;
	
	int refSpeciesIndex=0;
	
	int reflen;
	boolean ready;	

	KmerHasher columnHasher;

	public MultipleAlignmentInputSequence(String refSpecies, RootedBinaryPhylogeneticTree tree) {
		this.refSpecies = refSpecies;
		this.tree = tree;
		nSpecies = tree.getNumSpecies();
	}

	/** constructs a multiple alignment input sequence.
	 * @param speciesNames a list of the names of the species in the alignments
	 * @param consensuses a list of the consensus sequences for each species.  The consensus sequences shoudl form a multiple alignment, including gaps.
	 * @param tree a tree of all the species in the alignment with branch lengths.
	 */
	public MultipleAlignmentInputSequence(List<String> speciesNames, List<String> consensuses, String refSpecies, RootedBinaryPhylogeneticTree tree) {
		this(refSpecies, tree);
		setSpeciesAndConsensuses(speciesNames, consensuses);
	}

	public int getNumSpecies() {
		return nSpecies;
	}
	
	void setSpeciesAndConsensuses(List<String> speciesNames, List<String> consensuses) {
		this.speciesNames = speciesNames;
		this.consensuses = consensuses;

		// Determine number of species in multiple alignment;
		Assert.a(speciesNames.size() == consensuses.size());
		Assert.a(nSpecies == speciesNames.size());
		Assert.a(nSpecies >= 1,"Number of species was " + nSpecies + " and supposed to be >= 1");

		refSpeciesIndex = speciesNames.indexOf(refSpecies);
		Assert.a(refSpeciesIndex != -1, "Reference species ",refSpecies," not found in ",StringUtils.join(speciesNames.iterator(), ','));
		
		columnHasher = new KmerHasher(KmerHasher.ACGTother, nSpecies);
		 
		// Determine length of consensus (padded) sequence in multiple alignment
		consensusLength = consensuses.get(refSpeciesIndex).length();
		for (int spec = 0; spec<nSpecies; spec++) {
			Assert.a(consensuses.get(spec).length() == consensusLength);
		}
		
		h = new KmerHasher(KmerHasher.ACGTN,1);
		reflen = 0;
		ref2con = new ArrayList<Integer>();
		String refCon = consensuses.get(refSpeciesIndex);
		for (int cpos=0; cpos<consensusLength; cpos++) {
			if (h.hashable(refCon.charAt(cpos))) {
				ref2con.add(cpos);
				reflen++;
			}
		}
		Assert.a(ref2con.size() == reflen);

		
		con2refLeft = new int[consensusLength]; // at cpos, = max(0 , argmax_rpos( ref2con(rpos)< cpos ) ) 
		int refLeft = 0;
		for (int cpos=0; cpos<consensusLength; cpos++) {
			if (refLeft < (reflen-1)) {
				if (ref2con.get(refLeft+1) < cpos)  { refLeft++; }
			}
			con2refLeft[cpos] = refLeft;
		}
		con2refRight = new int[consensusLength]; // at cpos = min(reflen-1, argmin_rpos( ref2con(rpos)>cpos) )
		int refRight = reflen-1;
		for (int cpos=consensusLength-1; cpos>=0; cpos--) {
			if (refRight>0) {
				if (ref2con.get(refRight-1) > cpos) { refRight--; }
			}
			con2refRight[cpos] = refRight;
		}
		
		log.debug("consensus length = " + consensusLength + "    Reference length = " + reflen);
	}
	
	public MultipleAlignmentColumn getX(int ix) {
		MultipleAlignmentColumn ret = new MultipleAlignmentColumn(ix);
		return ret;
	}

	public int length() {
		return reflen;
	}

	public InputSequence<?> getComponent(String name) {
		throw new UnsupportedOperationException();
	}

	public Collection<String> listComponents() {
		throw new UnsupportedOperationException();
	}

	public int con2refLeft( int cpos ) {
		return con2refLeft[cpos];
	}

	public int con2refRight( int cpos ) {
		return con2refRight[cpos];
	}

	public int ref2con(int pos) {
		return ref2con.get(pos);
	}

	public char characterInPaddedAlignment(int consensusPosition, int speciesNumber) {
		return consensuses.get(speciesNumber).charAt(consensusPosition);
	}

	public int numSpecies() {
		return nSpecies;
	}

	public InputSequence<MultipleAlignmentColumn> subSequence(int start, int end) {
		Assert.a(start >= 1);
		Assert.a(end <= this.length());
		Assert.a(start <= end);
		
		int constart1 = ref2con(start-1)+1; // start coord in consensus, one-based inclusive
		int conend1 = ref2con(end-1)+1; // end coord in consensus, one-based inclusive
		
		ArrayList<String> newcon = new ArrayList<String>();
		
		for (int j=0; j<nSpecies; j++) {
			newcon.add(consensuses.get(j).substring(constart1-1,conend1));
		}
		
		MultipleAlignmentInputSequence MA = new MultipleAlignmentInputSequence(speciesNames,newcon,refSpecies, tree);
		
		return MA;
	}

	public List<String> getSpeciesNames()
	{
		return speciesNames;
	}
	
	public List<String> getConsensusSeqs()
	{
		return consensuses;
	}

	public int getColumnUniqueHash(int conpos) {
		int ret = 0;
		for(int i = 0; i<nSpecies; ++i) {
			ret = columnHasher.shiftHash(consensuses.get(i).charAt(conpos), ret);
		}
		return ret;
	}

	public int getConsensusLength() {
		return consensusLength;
	}

	public RootedBinaryPhylogeneticTree getTree() {
		return tree;
	}
	
	public PhylogeneticTreeFelsensteinOrder getFelsensteinOrder() {
		int n = speciesNames.size();
		String[] sn = new String[n];
		for (int j=0; j<n; j++) {
			sn[j] = speciesNames.get(j);
		}
		return tree.subtree(sn).getFelsensteinOrder(sn);
	}

	/** represents the column of the alignment at a given position on the reference sequence */
	public class MultipleAlignmentColumn {
		private int pos;
		private int cpos;

		/** constructs the column at this position */
		public MultipleAlignmentColumn(int pos) {
			this.pos = pos;
			cpos = ref2con(pos);
		}

		/** returns the multiple alignment input sequence that this column comes from.
		 * @return the owning input sequence
		 */
		public MultipleAlignmentInputSequence getMultipleAlignment() {
			return MultipleAlignmentInputSequence.this;
		}
		
		/** returns the number of species in the alignment
		 * @return the number of species in the alignment
		 */
		public int numSpecies() {
			return MultipleAlignmentInputSequence.this.numSpecies();
		}

		/** returns the value of this position in the alignment for the given species
		 * @param spec the species to retrieve.  The value will be the index of the species in the <code>speciesName</code> list.
		 * @return the character for this species in this column
		 */
		public char nucleotide(int spec) {
			return MultipleAlignmentInputSequence.this.characterInPaddedAlignment(cpos,spec);
		}

		/** returns a hash value for the characters in this column
		 * @return the hash value for the column
		 */
		public int getUniqueHash() {
			return MultipleAlignmentInputSequence.this.getColumnUniqueHash(pos);
		}

		/** returns the tree for this alignment
		 * @return the species tree for the alignment
		 */
		public RootedBinaryPhylogeneticTree getTree() {
			return MultipleAlignmentInputSequence.this.getTree();
		}
	}
}
