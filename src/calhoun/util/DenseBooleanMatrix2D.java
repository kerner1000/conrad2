/**
 * 
 */
package calhoun.util;

import java.io.Serializable;
import java.util.Arrays;

public class DenseBooleanMatrix2D implements Serializable {
	private static final long serialVersionUID = 8033567052720988542L;
	boolean[] array;
	int length;
	
	public DenseBooleanMatrix2D(int x, int y) {
		array = new boolean[x*y];
		length = y;
	}
	
	public boolean getQuick(int x, int y) {
		return array[x*length+y];
	}

	public void setQuick(int x, int y, boolean val) {
		array[x*length+y] = val;
	}
	
	public void assign(boolean val) {
		for(int i = 0; i<array.length; ++i) {
			array[i]= val;
		}
	}
	
	@Override
	public String toString() {
		StringBuffer b = new StringBuffer();
		b.append('\n');
		for(int i = 0; i<array.length/length;++i) {
			int row = i*length;
			for(int j = 0; j<length; ++j) {
				b.append(array[row+j]).append("   ");
			}
			b.append('\n');
		}
		return b.toString();
	}

	@Override
	public boolean equals(Object rhs) {
		DenseBooleanMatrix2D matrix = (DenseBooleanMatrix2D ) rhs;
		
		return Arrays.equals(array, matrix.array);
	}
}