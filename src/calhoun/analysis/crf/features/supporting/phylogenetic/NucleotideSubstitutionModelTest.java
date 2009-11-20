package calhoun.analysis.crf.features.supporting.phylogenetic;

import calhoun.util.ColtUtil;

public class NucleotideSubstitutionModelTest {

	public static void main(String[] args) {

		NucleotideSubstitutionModel M = new Kimura80Model(new double[]{1.2,0.5});
		
		System.out.println("Transition matrix at time 0.1:");
		System.out.println(ColtUtil.format(M.transitionMatrix(0.1)));
		
		System.out.println("Transition matrix at time 1.0:");
		System.out.println(ColtUtil.format(M.transitionMatrix(1.0)));
		
		System.out.println("Transition matrix at time 10:");
		System.out.println(ColtUtil.format(M.transitionMatrix(10)));
	}
}
