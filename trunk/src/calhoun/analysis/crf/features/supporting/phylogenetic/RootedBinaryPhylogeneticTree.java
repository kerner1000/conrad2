package calhoun.analysis.crf.features.supporting.phylogenetic;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.util.Assert;

public class RootedBinaryPhylogeneticTree implements Serializable {
	private static final long serialVersionUID = 4800435332663788163L;
	private static final Log log = LogFactory.getLog(RootedBinaryPhylogeneticTree.class);

	// Member variables
	public	ArrayList<BinaryTreeNode> T;
	
	// Constructors
	public RootedBinaryPhylogeneticTree( String ss ) {
		ss = ss.trim();
		weaklyValidateNewickString(ss);		
		T = new ArrayList<BinaryTreeNode>();
		growBranch(T,-1,ss.substring(0,ss.length()-1)); 
	}

	private RootedBinaryPhylogeneticTree( ArrayList<BinaryTreeNode> TT ) {      
		T = TT;
	}
	
	//////////////////////////////////////////////////
	// Public functions:
	public String newick() {
		Assert.a( T.get(0).p == -1 );
		String s = newick_recursion( 0 );
		s = s + ";";
		return s;
	}	

	public void summarize_tree() {
		log.debug("Now printing an extended summary of a RootedBiinaryPhylogeneticTree:");
		log.debug("Representation as a Newick String --->" + newick() );
		log.debug( "Total branch length (normalized): " + ( total_branch_length() ) );
		log.debug( "Longest branch length (normalized): " + ( longest_branch_length() ) );  
		for (int j=0; j<T.size(); j++) {
			log.debug( "" +  j + " -- " + T.get(j).toString() );
		}
	}
	
	public double total_branch_length() {
		double dist = 0;
		for (int j=0; j<T.size(); j++) {
			if (T.get(j).p == -1)
				Assert.a (T.get(j).d == 0);
			dist += T.get(j).d;
		}
		return dist;
	}
		
	public double longest_branch_length() {
		double longest = 0;
		for (int j=0; j<T.size(); j++) {
			if (T.get(j).p == -1) {
				double x = T.get(T.get(j).l).d + T.get(T.get(j).r).d;
				if (x>longest) longest = x;
			} else{
				double x = T.get(j).d;
				if (x>longest) longest = x;
			}
		}
		return longest;
	}

	public int getNumSpecies() {
		return getSpeciesSet().size();
	}
	
	/** Given an ordering of species in a multiple alignment file,
	 * determine an order for computations by which Felsenstein's algorithm can be performed, and the branch lengths involved at each step. 
	 * For example, for the phylogenetic tree</p>
	 * (((cnDT,cnDS),cnAB),(cnBB,cnBV))
	 * and the MSA in the order  (0=cnDT,1=cnDS,2=cnAB,3=cnBB,4=cnBV), can compute in this order:</p>
	 * 5 = combine(0,1)</p>
	 * 6 = combine(2,5)</p>
	 * 7 = combine(3,4)</p>
	 * 8 = combine(6,7)</p>
	 * and the final answer is at node 8.</p>
	 * 
	 * The computation for nodes 5-8 can be represented by</p>
	 * 		ileft = new int[]{0,5,3,6};</p>
	 *		iright = new int[]{1,2,4,7};</p>
	 * 
	 * @param msaOrder the (n) species presented in a multiple alignment, in order of their listing in the file
	 * @return  an ordering in which the (n-1) recursive calculations for Felsenstein's algorithm can be performed.
	 */
	public PhylogeneticTreeFelsensteinOrder getFelsensteinOrder(String[] msaOrder) {
		
		log.debug("Attempting to determine the order for a Felsenstein recursion given a phylogenetic tree");
		
		int nSpecies = msaOrder.length;
		int nSteps = nSpecies - 1;
				

		
		// Step 1: Make sure the species in msaOrder exactly correspond to the set of species in Phylogenetic Tree.
		Set<String> treeSpecies = getSpeciesSet();
		
		if (treeSpecies.size() != msaOrder.length) {
			Assert.a(false, "treeSpecies.size = " + treeSpecies.size() + "   and msaOrder.length = " + msaOrder.length);
		}
		for (int j=0; j<msaOrder.length; j++) {
			Assert.a(treeSpecies.contains(msaOrder[j]));
		}
		
		// Step 2: Initialize mappings from new to old and from old to new, such that for i=0..(nSpecies-1),
		//   T.get(new2old[i]).n = msaOrder[i]    and    old2new[new2old[i]] = i
		//   and also initilize the Boolean vector hasParent
		Map<Integer,Integer> new2old    = new HashMap<Integer,Integer>();
		Map<Integer,Integer> old2new    = new HashMap<Integer,Integer>();
		Map<Integer,Boolean> needParent = new HashMap<Integer,Boolean>();
		for (int j=0; j<nSpecies; j++) {
			int oldIndex = getSpeciesIndex(msaOrder[j]);
			log.debug("Species " + j + " is " + msaOrder[j]);
			old2new.put(oldIndex,j);
			new2old.put(j,oldIndex);
			needParent.put(j,true);
		}

		
		for (int step=0; step<nSteps; step++) {
			// Step 3: Find a nextNode to add on.  It must be:
			//   a) the parent of one of the nodes already in the new list which needs a parent
			//   b) both childern of nextNode are already in new list and need a parent
			int nextNodeOldIndex = -1;
			for (int j=0; j<new2old.size(); j++) {
				if (!needParent.get(j)) { continue; }
				int thisOldIndex = new2old.get(j);
				int parentOldIndex = T.get(new2old.get(j)).p;
				int leftOldIndex = T.get(parentOldIndex).l;
				int rightOldIndex = T.get(parentOldIndex).r;
				Assert.a( (thisOldIndex==leftOldIndex) || (thisOldIndex==rightOldIndex) ); // We went up and then came down two different ways, one of which must have resulted in no net movement
				if (!old2new.containsKey(leftOldIndex)) { continue; }
				if (!old2new.containsKey(rightOldIndex)) { continue; }	// Are both children of the considered parent already on the new list?
				nextNodeOldIndex = parentOldIndex;
				break;
			}
			Assert.a(nextNodeOldIndex != -1); // The process above has to have resulted in the identification of a suitable parent.
			
			// Step 4: Add this node to the mappings
			new2old.put(nSpecies+step,nextNodeOldIndex);
			old2new.put(nextNodeOldIndex,nSpecies+step);
			needParent.put(nSpecies+step,true);
			needParent.put(old2new.get(T.get(nextNodeOldIndex).l),false);
			needParent.put(old2new.get(T.get(nextNodeOldIndex).r),false);
		}
		
		// Step 5: Now that we know the mapping from new nodes to old nodes, just write down the order of computations.
		int[] ileft = new int[nSteps];
		int[] iright = new int[nSteps];
		double[] bleft = new double[nSteps];
		double[] bright = new double[nSteps];		
		
		for (int step=0; step <nSteps; step++) {
			int oldNodeIndex = new2old.get(nSpecies+step);
			ileft[step] = old2new.get(T.get(oldNodeIndex).l);
			iright[step] = old2new.get(T.get(oldNodeIndex).r);
			bleft[step] = T.get(new2old.get(ileft[step])).d;
			bright[step] = T.get(new2old.get(iright[step])).d;				
		}
		
//		for (int step=0; step <nSteps; step++) {
//			int oldNodeIndex = old2new.get(nSpecies+step);
//			ileft[step] = T.get(oldNodeIndex).l;
//			iright[step] = T.get(oldNodeIndex).r;
//			bleft[step] = T.get(old2new.get(ileft[step])).d;
//			bright[step] = T.get(old2new.get(iright[step])).d;				
//		}
		
		// Step 6: Construct the actual computation-order object, and return it.
		return new PhylogeneticTreeFelsensteinOrder( ileft, iright, bleft, bright);		
	}

	
	public RootedBinaryPhylogeneticTree subtree(String[] sn) {
		HashMap<String,Integer>  selected = new HashMap<String,Integer>();;
		
		for (int j=0; j<sn.length; j++) {
			selected.put(sn[j],0);
		}
		
		return subtree(selected);
	}
	
	
	private RootedBinaryPhylogeneticTree subtree( HashMap<String,Integer>  selected ) {
		ArrayList<BinaryTreeNode> oldT = T;
		
		for (int j=0; j<oldT.size(); j++) {
			if (oldT.get(j).l == -1) {
				Assert.a( oldT.get(j).r == -1 );
				Integer a = selected.get(oldT.get(j).n);
				if (a != null) {
					selected.put(oldT.get(j).n ,   selected.get(oldT.get(j).n).intValue() +1  );
					oldT.get(j).lm=true;
					oldT.get(j).rm=true;
					int s = j;  int p = oldT.get(s).p;
					while (p != -1) {
						Assert.a ( (oldT.get(p).l == s) || (oldT.get(p).r == s) );
						if (oldT.get(p).l == s) { oldT.get(p).lm = true; }
						if (oldT.get(p).r == s) { oldT.get(p).rm = true; }
						s = p; p = oldT.get(s).p;
					}
				}
			}
		}
		
		log.debug( "About to doublecheck that all selected species were found exactly once...");
		Iterator ii = selected.entrySet().iterator();
		
		while ( ii.hasNext() ) {
			Map.Entry<String,Integer> me = (Entry<String, Integer>) ii.next();
			if (me.getValue() != 1) {
				log.debug("  NOT FOUND UNIQUELY: " + me.getKey() +  "   " + me.getValue() );
				
			}
			
		}
		log.debug( "DONE" );		
		
		
		
		Map<Integer,Integer> old_to_new = new HashMap<Integer,Integer>();
		int numNew = 0;
		for (int j=0; j<oldT.size(); j++ ) {
			if (oldT.get(j).lm && oldT.get(j).rm) {
				old_to_new.put(j,numNew);
				numNew++;
			}
		}
		
		ArrayList<BinaryTreeNode> newT = new ArrayList<BinaryTreeNode>();
		for (int j=0; j<numNew; j++) {
			newT.add(new BinaryTreeNode());
		}
		
		for (int j=0; j<oldT.size(); j++) {
			if (!oldT.get(j).lm || !oldT.get(j).rm) continue;
			int s = old_to_new.get(j);
			newT.set(s , oldT.get(j) );
			
			int p = oldT.get(j).p;  double dist=oldT.get(j).d;
			while ( p!= -1 ) {
				if (oldT.get(p).lm && oldT.get(p).rm) break;
				dist += oldT.get(p).d;
				p = oldT.get(p).p;
			}
			if (p==-1) {
				newT.get(s).p=-1; newT.get(s).d=0;
			} else {
				newT.get(s).p=old_to_new.get(p);  newT.get(s).d=dist;
			}
			
			
			int l = oldT.get(j).l;
			if (l == -1) {
				newT.get(s).l = -1;
			} else {
				while (!oldT.get(l).lm || !oldT.get(l).rm) {
					Assert.a (oldT.get(l).lm || oldT.get(l).rm);
					if (oldT.get(l).lm) { l = oldT.get(l).l; continue; }
					if (oldT.get(l).rm) { l = oldT.get(l).r; continue; }
				}
				newT.get(s).l = old_to_new.get(l);
			}
			
			int r = oldT.get(j).r;
			if (r == -1) {
				newT.get(s).r = -1;
			} else {
				while (!oldT.get(r).lm || !oldT.get(r).rm) {
					Assert.a (oldT.get(r).lm || oldT.get(r).rm);
					if (oldT.get(r).lm) { r = oldT.get(r).l; continue; }
					if (oldT.get(r).rm) { r = oldT.get(r).r; continue; }
				}
				newT.get(s).r = old_to_new.get(r);
			}
		}
		
		RootedBinaryPhylogeneticTree subRBPT = new RootedBinaryPhylogeneticTree(newT);
		return subRBPT;
	}
	
	/////////////////////////////////////////////////////////
	// Internal private functions:

//	private Boolean containsSpecies(String name) {
//		Boolean ret = false;
//		for (int j=0; j<T.size(); j++) {
//			if (T.get(j).n == name) {
//				ret = true;
//			}
//		}
//		return ret;
//	}

	private Integer getSpeciesIndex(String name) {
		Integer ret = -1;
		for (int j=0; j<T.size(); j++) {
			String temp =  T.get(j).n;
			//System.out.println("This name is " + temp);
			if ( name.equals(temp) ) {
				ret = j;
				//System.out.println("  Species " + j + " = " + name);
			} else {
				//System.out.println("  --->" + name + "<--- and --->" + temp + "<--- are not equal");
			}
		}
		Assert.a(ret != -1," Could not find species " + name);
		return ret;
	}
	
	public Set<String> getSpeciesSet() {
		Set<String> ret = new HashSet<String>();
		for (int j=0; j<T.size(); j++) {
			BinaryTreeNode btn = T.get(j);
			if (btn.l == -1) { // this node has no left child
				Assert.a(btn.r == -1); // if no left child, then its a leaf and hence also no right child either.
				String name = btn.n;
				Assert.a(name != ""); // I require that each leaf node have a name.
				Assert.a(!ret.contains(name)); // I want no repitition; if same species is listed in tree twice then fail here
				ret.add(name);
			}
		}
		return ret;
	}
	
	private String newick_recursion( int top ) {
		String s = "";  
		
		if (T.get(top).l == -1) {
			Assert.a( T.get(top).r == -1);
			s = s + T.get(top).n;
		} else {
			s = s + "(";
			s = s + newick_recursion( T.get(top).l );
			s = s + ",";
			s = s + newick_recursion( T.get(top).r );
			s = s + ")";
		}
		s = s + ":" + T.get(top).d;
		return s;
	};
	
	private void growBranch( ArrayList<BinaryTreeNode> T1, int parent, String S ) {
		
		int currentNode = T1.size();
		
		int x = S.lastIndexOf(":");
		
		Assert.a( x != -1, "we require that all nodes in Newick format have a branch length" );  
		Assert.a( x != S.length()-1 );
		
		double dist = Double.parseDouble(S.substring(x+1));
		
		String S2 = S.substring(0,x);
		String S3;
		
		if (S2.charAt(0)=='(') {
			Assert.a (S2.charAt(S2.length()-1) == ')');
			S3 = S2.substring(1,S2.length()-1);
		} else {
			S3 = S2;
		}
		
		int depth=0;
		for (int j=0; j<(S3.length()-1); j++) {
			if (S3.charAt(j) == '(') depth++;
			if (S3.charAt(j) == ')') depth--;
			if ((depth == 0) && (S3.charAt(j)==',')) {
				String leftString = S3.substring(0,j);
				String rightString = S3.substring(j+1,S3.length());
				BinaryTreeNode TN = new BinaryTreeNode( parent, -1,-1, dist, "");
				T1.add(TN);
				T1.get(currentNode).l = T1.size();
				growBranch(T1,currentNode,leftString);
				T1.get(currentNode).r = T1.size();
				growBranch(T1,currentNode,rightString);
				return;
			}
		}
		Assert.a(depth==0);
		
		BinaryTreeNode T2 = new BinaryTreeNode( parent,-1,-1,dist,S3);
		T1.add(T2);
		return;
	}	
	
	private void weaklyValidateNewickString(String ss) {
		int depth=0;
		Assert.a( ss.charAt(ss.length()-1) == ';', "Invalid newick tree: "+ss);
		for (int j=0; j<(ss.length()-1); j++) {
			if (ss.charAt(j) == '(') depth++;
			if (ss.charAt(j) == ')') depth--;
			Assert.a(depth >= 0);
		}
		Assert.a(depth==0);
	}

	public int nSpecies() {
		int m = T.size(); // should be m=2*n-1, where n is number of species
		Assert.a((m%2)==1);
		return m/2;
	}

	
}
