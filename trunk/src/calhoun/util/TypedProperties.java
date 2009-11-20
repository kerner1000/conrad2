package calhoun.util;

import java.util.Properties;

import org.apache.commons.lang.StringUtils;

public class TypedProperties extends Properties {
	private static final long serialVersionUID = 8275327757663063251L;

	public TypedProperties() {
	}
	
	public TypedProperties(Properties props) {
		super(props);
	}
	
	public int getIntProperty(String name, int defValue) {
		String val = getProperty(name);
		return (val == null) ? defValue : Integer.parseInt(val);
	}

	public double getDoubleProperty(String name, double defValue) {
		String val = getProperty(name);
		return (val == null) ? defValue : Double.parseDouble(val);
	}

	public double[] getDoubleArrayProperty(String name, double[] defValue) {
		String val = getProperty(name);
		if(val == null) {
			return defValue;
		}
		String[] vals = StringUtils.split(val, ',');
		double[] ret = new double[vals.length];
		for(int i = 0; i<vals.length; ++i) {
			ret[i] = Double.parseDouble(vals[i]);
		}
		return ret;
	}

	public int[] getIntArrayProperty(String name, int[] defValue) {
		String val = getProperty(name);
		if(val == null) {
			return defValue;
		}
		String[] vals = val.split(",");
		int[] ret = new int[vals.length];
		for(int i = 0; i<vals.length; ++i) {
			ret[i] = Integer.parseInt(vals[i]);
		}
		return ret;
	}

	public short[] getShortArrayProperty(String name, short[] defValue) {
		String val = getProperty(name);
		if(val == null) {
			return defValue;
		}
		String[] vals = val.split(",");
		short[] ret = new short[vals.length];
		for(int i = 0; i<vals.length; ++i) {
			ret[i] = Short.parseShort(vals[i]);
		}
		return ret;
	}

	public boolean  getBooleanProperty(String name, boolean defValue) {
		String val = getProperty(name);
		return (val == null) ? defValue : Boolean.parseBoolean(val);
	}

}
