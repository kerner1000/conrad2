package calhoun.analysis.crf.executables.test;

public class TestPreincrement {

	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
	
		System.out.println("-4 mode 3 is " + ((-4 %3+9)%3));
		
		System.out.println("Maximum value of a short is " + Short.MAX_VALUE);
		System.out.println("Minimum value of a short is " + Short.MIN_VALUE);

		System.out.println("Maximum value of a Integer is " + Integer.MAX_VALUE);
		System.out.println("Minimum value of a Integer is " + Integer.MIN_VALUE);
		
		int j;
		for (j=0; j<4; j++) {
			System.out.println(j);
		}
		System.out.println("  "+j);
		for (j=0; j<4; ++j) {
			System.out.println(j);
		}
		System.out.println("  "+j);
		for (j=3; j>0; j--) {
			System.out.println(j);
		}
		System.out.println("  "+j);
		for (j=3; j>0; --j) {
			System.out.println(j);
		}
		System.out.println("  "+j);
		
		int type = 2;
		System.out.println(" type==2 ? 3:1 -----> " + (type == 2 ? 3:1) );
		
		
		String s = "(164716,215910,-,1.0)(194048,218199,-,1.0)";	
		
		String[] fields = s.split("[]");
		for (int k=0; k<fields.length; k++) {
			System.out.println("-->" + fields[k] + "<--");
		}
	
		String[] rings = fields[1].split("[,)]");
		for (int k=0; k<rings.length; k++) {
			System.out.println("---->" + rings[k] + "<--");
		}
		
		
	}

}
