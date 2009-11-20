package calhoun.analysis.crf.features.interval13;

import java.util.Arrays;

import calhoun.analysis.crf.BeanModel;


public class Interval13Model extends BeanModel {
	private static final long serialVersionUID = 3959312826759045449L;

	private static boolean narrowBoundaries = false;
	
	public void setNarrowBoundaries(boolean narrowBoundaries) {
		Interval13Model.narrowBoundaries = narrowBoundaries;
	}
	
	public static int getPadExon5prime() {
		if (narrowBoundaries) {
			return 4;
		} else {
			return 6;
		}
	}

	public static int getPadExon3prime() {
		if (narrowBoundaries) {
			return 1;
		} else {
			return 3;
		}
	}
	
	public static int getPadIntron5prime() {
		return 6;
	}	
	
	public static int getPadIntron3prime() {
		return 9;
	}	
	
	public static int getPadIntergenic() {
		return 9;
	}			
	
	
	static Node[] nodeArray = new Node[] {
		new Node(0, "intergenic"),
		new Node(1, "exon0"),
		new Node(2, "exon1"),
		new Node(3, "exon2"),
		new Node(4, "intron0"),
		new Node(5, "intron1"),
		new Node(6, "intron2"),
		new Node(7, "exon0m"),
		new Node(8, "exon1m"),
		new Node(9, "exon2m"),
		new Node(10, "intron0m"),
		new Node(11, "intron1m"),
		new Node(12, "intron2m"),
	};
	
	static Edge[] edgeArray = new Edge[] {
		new Edge(nodeArray[0], nodeArray[0]),
		new Edge(nodeArray[1], nodeArray[1]),
		new Edge(nodeArray[2], nodeArray[2]),
		new Edge(nodeArray[3], nodeArray[3]),
		new Edge(nodeArray[4], nodeArray[4]),
		new Edge(nodeArray[5], nodeArray[5]),
		new Edge(nodeArray[6], nodeArray[6]),
		new Edge(nodeArray[7], nodeArray[7]),
		new Edge(nodeArray[8], nodeArray[8]),
		new Edge(nodeArray[9], nodeArray[9]),
		new Edge(nodeArray[10], nodeArray[10]),
		new Edge(nodeArray[11], nodeArray[11]),
		new Edge(nodeArray[12], nodeArray[12]),
		new Edge(nodeArray[0], nodeArray[1]),
		new Edge(nodeArray[0], nodeArray[2]),
		new Edge(nodeArray[0], nodeArray[3]),
		new Edge(nodeArray[0], nodeArray[7]),
		new Edge(nodeArray[0], nodeArray[8]),
		new Edge(nodeArray[0], nodeArray[9]),
		new Edge(nodeArray[1], nodeArray[0]),
		new Edge(nodeArray[2], nodeArray[0]),
		new Edge(nodeArray[3], nodeArray[0]),
		new Edge(nodeArray[7], nodeArray[0]),
		new Edge(nodeArray[8], nodeArray[0]),
		new Edge(nodeArray[9], nodeArray[0]),
		new Edge(nodeArray[1], nodeArray[4]),
		new Edge(nodeArray[1], nodeArray[5]),
		new Edge(nodeArray[1], nodeArray[6]),
		new Edge(nodeArray[2], nodeArray[4]),
		new Edge(nodeArray[2], nodeArray[5]),
		new Edge(nodeArray[2], nodeArray[6]),
		new Edge(nodeArray[3], nodeArray[4]),
		new Edge(nodeArray[3], nodeArray[5]),
		new Edge(nodeArray[3], nodeArray[6]),
		new Edge(nodeArray[4], nodeArray[1]),
		new Edge(nodeArray[4], nodeArray[2]),
		new Edge(nodeArray[4], nodeArray[3]),
		new Edge(nodeArray[5], nodeArray[1]),
		new Edge(nodeArray[5], nodeArray[2]),
		new Edge(nodeArray[5], nodeArray[3]),
		new Edge(nodeArray[6], nodeArray[1]),
		new Edge(nodeArray[6], nodeArray[2]),
		new Edge(nodeArray[6], nodeArray[3]),
		new Edge(nodeArray[7], nodeArray[10]),
		new Edge(nodeArray[7], nodeArray[11]),
		new Edge(nodeArray[7], nodeArray[12]),
		new Edge(nodeArray[8], nodeArray[10]),
		new Edge(nodeArray[8], nodeArray[11]),
		new Edge(nodeArray[8], nodeArray[12]),
		new Edge(nodeArray[9], nodeArray[10]),
		new Edge(nodeArray[9], nodeArray[11]),
		new Edge(nodeArray[9], nodeArray[12]),
		new Edge(nodeArray[10], nodeArray[7]),
		new Edge(nodeArray[10], nodeArray[8]),
		new Edge(nodeArray[10], nodeArray[9]),
		new Edge(nodeArray[11], nodeArray[7]),
		new Edge(nodeArray[11], nodeArray[8]),
		new Edge(nodeArray[11], nodeArray[9]),
		new Edge(nodeArray[12], nodeArray[7]),
		new Edge(nodeArray[12], nodeArray[8]),
		new Edge(nodeArray[12], nodeArray[9]),
	};
	
	public Interval13Model() {
		nodes = Arrays.asList(nodeArray);
		edges = Arrays.asList(edgeArray);
	}
}
