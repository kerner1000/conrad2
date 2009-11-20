package calhoun.analysis.crf.executables.test;

import calhoun.util.Assert;

public class SwitchStatementTest {

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {	

		int x = 3;
	
		
		
		System.out.println("Set x = " + x);
		
		switch (x) {
		case 1: 
			System.out.println("In this case, x = 1"); 
			break;
		case 2:
			System.out.println("In this case, x = 2"); 
			break;
		case 3:
			System.out.println("In this case, x = 3"); 
			break;
		default:
			Assert.a(false,"Unanticipated case for x");	
		}

		int[] a = new int[]{3,47,1,3};
		System.out.println(" before, a[0] = " + a[0]);
		doesThisModify(a);
		System.out.println(" after, a[0] = " + a[0]);		
	
		for (int y : a) {
			System.out.println(" after, a[?] = " + y);				
			
		}
		
		return;
	}

	public static void doesThisModify(int[] a) {
		a[0]++;
	}
	
}
