/*
 * Created on Nov 19, 2003
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package calhoun.util;

import java.beans.PropertyDescriptor;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import junit.framework.TestCase;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

/** Basic test case with file matching (XML and regular) added as well as improved logging for errors. 
 */
public abstract class AbstractTestCase extends TestCase {
	private static final Log log = LogFactory.getLog(AbstractTestCase.class);

	public AbstractTestCase() {
    }
	
	public AbstractTestCase(String arg0) {
		super(arg0);
	}

	/** Overrides the default runTest to provide error logging.  Using this runBare logs setUp and tearDown errors.
	 * 
	 */
	@Override
	public void runBare() throws Throwable {
		try {
			log.info("Running Test: " + this.getClass().getName() + " " + getName());
			setUp();
		} catch (Throwable ex) {
			log.error(getName() + ": ", ex);
			throw ex;
		}
		try {
			runTest();
		} catch (Throwable ex) {
			log.error(getName() + ": ", ex);
			throw ex;
		} finally {
			tearDown();
		}
	}

	
	private void failCompare(String name, Object expect, Object actual) throws Exception {
		fail("Property compare failed on " + name + ".  Expected: " + expect + " Actual: " + actual);
	}

	private static Set basicTypes = new HashSet();
	static {
		basicTypes.add(String.class);
		basicTypes.add(Byte.class);
		basicTypes.add(Byte.TYPE);
		basicTypes.add(Boolean.class);
		basicTypes.add(Boolean.TYPE);
		basicTypes.add(Short.class);
		basicTypes.add(Short.TYPE);
		basicTypes.add(Integer.class);
		basicTypes.add(Integer.TYPE);
		basicTypes.add(Long.class);
		basicTypes.add(Long.TYPE);
		basicTypes.add(Character.class);
		basicTypes.add(Character.TYPE);
		basicTypes.add(Double.class);
		basicTypes.add(Double.TYPE);
		basicTypes.add(Float.class);
		basicTypes.add(Float.TYPE); 
	}

	public void assertXmlFilesMatch(String template, String actual) throws Exception {
		String result = XmlDiff.compareFiles(template, actual);
		log.debug("Result: "+result);
		if (result != null) {
			log.warn(result);
			fail();
		}
	}

	public void assertArrayEquals(byte[] a, byte[] b) {
		for(int i=0;i<a.length;i++) {
			assertEquals(a[i], b[i]);
		}
	}
	
	public void assertArrayEquals(double[] a, double[] b) {
		for(int i=0;i<a.length;i++) {
			assertEquals(a[i], b[i]);
		}
	}
	
	public void assertArrayEquals(double[] a, double[] b, double tolerance) {
		for(int i=0;i<a.length;i++) {
			assertEquals(a[i], b[i], tolerance);
		}
	}
	
	public void assertArrayEquals(int[] a, int[] b) {
		for(int i=0;i<a.length;i++) {
			assertEquals(a[i], b[i]);
		}
	}
	
	public void assertBeginsWith(byte[] a, byte[] b) {
		for(int i=0;i<a.length;i++) {
			assertEquals(a[i], b[i]);
		}
	}
	
	public void assertEndsWith(byte[] a, byte[] b) {
		for(int i=0;i<a.length;i++) {
			assertEquals(a[i], b[b.length-a.length+i]);
		}
	}
	
	public void assertEquals(byte[] a, byte[] b) {
		assertEquals(a.length, b.length);
		for(int i=0;i<a.length;i++) {
			assertEquals(a[i], b[i]);
		}
	}
	
	public void assertStringEquals(String [] expected, List l) {
		assertStringEquals(expected, (String[])l.toArray(new String[l.size()]));
	}
	
	public void assertStringEquals(String [] expected, String[] strings) {
		for (int i = 0; i < strings.length; ++i) {
			if (expected.length > i) {
				assertEquals(expected[i], strings[i]);
			}
		}
		assertEquals(expected.length, strings.length);
	}
	
	public void assertFilesMatch(String expected, String actual) throws Exception {
		FileDiff diff = new FileDiff(expected, actual);
		boolean mismatch = diff.execute();
		assertFalse(diff.toString(), mismatch);
	}

	public void assertFlatFileRecordsMatch(String expected, String actual) throws Exception {
		String s = FileUtil.fileRecordCompare(expected, actual);
		assertNull(s, s);
	}
}
