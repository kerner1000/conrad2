package calhoun.analysis.crf.features.supporting.phylogenetic;

import java.io.Serializable;

public class BinaryTreeNode implements Serializable{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -2109391655262800568L;
	public  int    p;   // parent;
	public  int    l;   // left child
	public  int    r;   // right child
	public  double d;   // branch distance
	public  String n;   // name of node
	public  Boolean   lm;  // left mark
	public  Boolean   rm;  // right mark
	
	
	BinaryTreeNode() {
		p  = -1;
		l  = -1;
		r  = -1;
		d  = 0.0;
		n  = "";
		lm = false;
		rm = false;
	}
	
	
	BinaryTreeNode(int parent, int left, int right, double dist, String name) {
		p  = parent;
		l  = left;
		r  = right;
		d  = dist;
		n  = name;
		lm = false;
		rm = false;
	}
	
	@Override
	public String toString() {
		String ret = "parent: " + p + "  left: " + l + "  right: " + r + "  distance: " + d + "  name: " + n;
		return ret;
	}
	
	
	
}
