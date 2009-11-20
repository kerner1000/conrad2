package calhoun.analysis.crf.features.interval29;

import java.util.Arrays;
import java.util.List;

import calhoun.analysis.crf.BeanModel;


public class Interval29Model extends BeanModel {
	private static final long serialVersionUID = 3959312826759045449L;

	private static boolean narrowBoundaries = false;
	
	public void setNarrowBoundaries(boolean narrowBoundaries) {
		Interval29Model.narrowBoundaries = narrowBoundaries;
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
		new Node(13, "ig-e"),
		new Node(14, "e-ig"),
		new Node(15, "e-i0"),
		new Node(16, "e-i1"),
		new Node(17, "e-i2"),
		new Node(18, "i-e0"),
		new Node(19, "i-e1"),
		new Node(20, "i-e2"),
		new Node(21, "ig-em"),
		new Node(22, "em-ig"),
		new Node(23, "em-i0m"),
		new Node(24, "em-i1m"),
		new Node(25, "em-i2m"),
		new Node(26, "im-e0m"),
		new Node(27, "im-e1m"),
		new Node(28, "im-e2m")
	};
	
	static Edge[] edgeArray = new Edge[] {
		// To self
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
		new Edge(nodeArray[13], nodeArray[13]),
		new Edge(nodeArray[14], nodeArray[14]),
		new Edge(nodeArray[15], nodeArray[15]),
		new Edge(nodeArray[16], nodeArray[16]),
		new Edge(nodeArray[17], nodeArray[17]),
		new Edge(nodeArray[18], nodeArray[18]),
		new Edge(nodeArray[19], nodeArray[19]),
		new Edge(nodeArray[20], nodeArray[20]),
		new Edge(nodeArray[21], nodeArray[21]),
		new Edge(nodeArray[22], nodeArray[22]),
		new Edge(nodeArray[23], nodeArray[23]),
		new Edge(nodeArray[24], nodeArray[24]),
		new Edge(nodeArray[25], nodeArray[25]),
		new Edge(nodeArray[26], nodeArray[26]),
		new Edge(nodeArray[27], nodeArray[27]),
		new Edge(nodeArray[28], nodeArray[28]),
		
		// Intergenic to intergenic-exon_i boundary
		new Edge(nodeArray[0], nodeArray[13]),
		// Intergenic-exon_i boundary to exon_i
		new Edge(nodeArray[13], nodeArray[1]),
		new Edge(nodeArray[13], nodeArray[2]),
		new Edge(nodeArray[13], nodeArray[3]),
		// Intergenic to intergenic-exon_im boundary
		new Edge(nodeArray[0], nodeArray[21]),
		// Intergenic-exonm boundary to exonm
		new Edge(nodeArray[21], nodeArray[7]),
		new Edge(nodeArray[21], nodeArray[8]),
		new Edge(nodeArray[21], nodeArray[9]),
		
		// Exon_i to exon-intergenic boundary
		new Edge(nodeArray[1], nodeArray[14]),
		new Edge(nodeArray[2], nodeArray[14]),
		new Edge(nodeArray[3], nodeArray[14]),
		// Exon-intergenic boundary to intergenic
		new Edge(nodeArray[14], nodeArray[0]),
		// Exon_im to exonm-intergenic boundary
		new Edge(nodeArray[7], nodeArray[22]),
		new Edge(nodeArray[8], nodeArray[22]),
		new Edge(nodeArray[9], nodeArray[22]),
		// Exonm-intergenic boundary to intergenic
		new Edge(nodeArray[22], nodeArray[0]),
		
		// Exon_i to exon-intron_j boundary
		new Edge(nodeArray[1], nodeArray[15]),
		new Edge(nodeArray[2], nodeArray[15]),
		new Edge(nodeArray[3], nodeArray[15]),
		new Edge(nodeArray[1], nodeArray[16]),
		new Edge(nodeArray[2], nodeArray[16]),
		new Edge(nodeArray[3], nodeArray[16]),
		new Edge(nodeArray[1], nodeArray[17]),
		new Edge(nodeArray[2], nodeArray[17]),
		new Edge(nodeArray[3], nodeArray[17]),
		// Exon-intron_i boundary to intron_i
		new Edge(nodeArray[15], nodeArray[4]),
		new Edge(nodeArray[16], nodeArray[5]),
		new Edge(nodeArray[17], nodeArray[6]),
		
		// Intron_i to intron-exon_j boundary
		new Edge(nodeArray[4], nodeArray[18]),
		new Edge(nodeArray[5], nodeArray[18]),
		new Edge(nodeArray[6], nodeArray[18]),
		new Edge(nodeArray[4], nodeArray[19]),
		new Edge(nodeArray[5], nodeArray[19]),
		new Edge(nodeArray[6], nodeArray[19]),
		new Edge(nodeArray[4], nodeArray[20]),
		new Edge(nodeArray[5], nodeArray[20]),
		new Edge(nodeArray[6], nodeArray[20]),	
		// Intron-exon_i boundary to exon_i
		new Edge(nodeArray[18], nodeArray[1]),
		new Edge(nodeArray[19], nodeArray[2]),
		new Edge(nodeArray[20], nodeArray[3]),

		// Exon_im to exonm-intron_jm boundary
		new Edge(nodeArray[7], nodeArray[23]),
		new Edge(nodeArray[8], nodeArray[23]),
		new Edge(nodeArray[9], nodeArray[23]),
		new Edge(nodeArray[7], nodeArray[24]),
		new Edge(nodeArray[8], nodeArray[24]),
		new Edge(nodeArray[9], nodeArray[24]),
		new Edge(nodeArray[7], nodeArray[25]),
		new Edge(nodeArray[8], nodeArray[25]),
		new Edge(nodeArray[9], nodeArray[25]),
		// Exonm-intron_im boundary to intron_im
		new Edge(nodeArray[23], nodeArray[10]),
		new Edge(nodeArray[24], nodeArray[11]),
		new Edge(nodeArray[25], nodeArray[12]),

		// Intron_im to intronm-exon_jm boundary
		new Edge(nodeArray[10], nodeArray[26]),
		new Edge(nodeArray[11], nodeArray[26]),
		new Edge(nodeArray[12], nodeArray[26]),
		new Edge(nodeArray[10], nodeArray[27]),
		new Edge(nodeArray[11], nodeArray[27]),
		new Edge(nodeArray[12], nodeArray[27]),
		new Edge(nodeArray[10], nodeArray[28]),
		new Edge(nodeArray[11], nodeArray[28]),
		new Edge(nodeArray[12], nodeArray[28]),
		// Intronm-exonm boundary to exonm
		new Edge(nodeArray[26], nodeArray[7]),
		new Edge(nodeArray[27], nodeArray[8]),
		new Edge(nodeArray[28], nodeArray[9])  // done
	};
	
	public List checkValidTransitions(List data) {
		return Interval29Tools.checkValidTransitions(data);
	}
	
	public Interval29Model() {
		nodes = Arrays.asList(nodeArray);
		edges = Arrays.asList(edgeArray);
	}
}
