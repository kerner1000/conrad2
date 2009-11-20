package calhoun.analysis.crf.test;

import java.util.Arrays;

import calhoun.analysis.crf.BeanModel;


public class ZeroOrderModel extends BeanModel { 
	private static final long serialVersionUID = 3959312826759045449L;

	static Node[] nodeArray = new Node[] {
			new Node(0, "lowGC"),
			new Node(1, "highGC")
	};
	
	static Edge[] edgeArray = new Edge[] {
			new Edge(nodeArray[0], nodeArray[0]),
			new Edge(nodeArray[0], nodeArray[1]),
			new Edge(nodeArray[1], nodeArray[1]),
			new Edge(nodeArray[1], nodeArray[0])
	};
	
	public ZeroOrderModel() {
		nodes = Arrays.asList(nodeArray);
		edges = Arrays.asList(edgeArray);
	}
}
