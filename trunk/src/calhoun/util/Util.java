/* The Broad Institute SOFTWARE COPYRIGHT NOTICE AGREEMENT
This software and its documentation are copyright 2003 by the Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
This software is supplied without any warranty or guaranteed support whatsoever. Neither the Broad Institute nor MIT can be responsible for its use, misuse, or functionality. 
*/
package calhoun.util;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

public class Util {

	/** Static helpers.  Do not instantiate */
	private Util() {
		super();
	}

	public static<T> void addAll(Collection<T> coll, Iterator<? extends T> iter) {
		while(iter.hasNext()) {
			coll.add(iter.next());
		}
	}
	
	public static void normalizeWeights(float[] weights) {
		float total = 0.0f;
		for(int i=0; i<weights.length; ++i) {
			total += weights[i];
		}
		for(int i=0; i<weights.length; ++i) {
			weights[i] = weights[i]/total;
		}
	}
	
	public static int[] convertIntList(List<Integer> list) {
		int[] labelArray = new int[list.size()];
		for(int j=0; j<labelArray.length; ++j) {
			labelArray[j] = list.get(j);
		}
		return labelArray;
	}

	public static float[] convertFloatList(List<Float> list) {
		float[] labelArray = new float[list.size()];
		for(int j=0; j<labelArray.length; ++j) {
			labelArray[j] = list.get(j);
		}
		return labelArray;
	}

	public static List<Double> convertDoubleArray(double[] list) {
		List<Double> labelArray = new ArrayList<Double>(list.length);
		for(int j=0; j<list.length; ++j) {
			labelArray.add(list[j]);
		}
		return labelArray;
	}

	static public boolean safeEquals(Object x, Object y) {
		if(x == null)
			return y == null;
		else
			return x.equals(y);		
	}

	/** make a copy of any serializable object graph */
	static public Object deepClone(Object o) {
		try { 
			Object copy;
			
			ByteArrayOutputStream baos = new ByteArrayOutputStream(); 
			ObjectOutputStream oos = new ObjectOutputStream(baos);
			oos.writeObject(o);
			oos.close();
			ByteArrayInputStream bios = new ByteArrayInputStream(baos.toByteArray());
			ObjectInputStream ois = new ObjectInputStream(bios);
			copy = ois.readObject();
			
			return copy;
		} catch (Exception ex) {
			throw new ErrorException("Could not copy object "+o.toString(), ex);
		}
	}
}
