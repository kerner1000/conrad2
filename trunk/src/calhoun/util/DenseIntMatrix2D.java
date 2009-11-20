/**
 * 
 */
package calhoun.util;

import java.io.Serializable;
import java.util.Arrays;

public final class DenseIntMatrix2D implements Serializable {
	private static final long serialVersionUID = 993972669435690618L;
	int[] array;
	int length;
	
	public DenseIntMatrix2D(int x, int y) {
		array = new int[x*y];
		length = y;
	}
	
	public final int getQuick(final int x, final int y) {
		return array[x*length+y];
	}

	public final void setQuick(final int x, final int y, final int val) {
		array[x*length+y] = val;
	}
	
	public final void assign(final int val) {
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
		DenseIntMatrix2D matrix = (DenseIntMatrix2D) rhs;
		
		return Arrays.equals(array, matrix.array);
	}
}
