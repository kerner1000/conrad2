package calhoun.analysis.crf.test;

import java.io.File;

import calhoun.analysis.crf.io.ExtensionMapper;
import calhoun.util.AbstractTestCase;

public class FileInputHandlerTest extends AbstractTestCase {
	
	public void testFilenameMapper() throws Exception {
		ExtensionMapper m = new ExtensionMapper();
		m.setExtension(".new");
		assertEquals(new File("test.new"), m.mapFilename(new File("test.old")));
		assertEquals(new File("/absolute/path/test.new"), m.mapFilename(new File("/absolute/path/test.old")));

		m.setAppend(true);
		assertEquals(new File("/absolute/path/test.old.new"), m.mapFilename(new File("/absolute/path/test.old")));
	}
}
