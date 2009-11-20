package calhoun.analysis.crf.test;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.Conrad;
import calhoun.util.AbstractTestCase;

public class RealIOTest extends AbstractTestCase {
	private static final Log log = LogFactory.getLog(RealIOTest.class);

	public void testRealFormatsNeg() throws Exception {
		Conrad.main(new String[] {"train", "test/input/realFormats/neg/config.xml", "test/input/realFormats/neg", "test/working/realFormats.ser"});
	}

	public void testRealFormatsMulti() throws Exception {
		Conrad.main(new String[] {"train", "test/input/realFormats/multiConfig.xml", "test/input/realFormats/aln", "test/working/realFormats.ser"});
	}

	public void testRealFormats2() throws Exception {
		Conrad.main(new String[] {"train", "test/input/realFormats/config2.xml", "test/input/realFormats", "test/working/realFormats.ser"});
		Conrad.main(new String[] {"test", "test/working/realFormats.ser", "test/input/realFormats", "test/working/realFormatsPredicted.txt"});
	}

	public void testRealFormats() throws Exception {
		Conrad.main(new String[] {"train", "test/input/realFormats/config.xml", "test/input/realFormats", "test/working/realFormats.ser"});
		Conrad.main(new String[] {"test", "test/working/realFormats.ser", "test/input/realFormats", "test/working/realFormatsPredicted.txt"});
	}

	public void testModulo() {
		log.warn("-1%3="+(-1%3));
		log.warn("-2%3="+(-2%3));
		log.warn("-3%3="+(-3%3));
		log.warn("-4%3="+(-4%3));
	}
	
}
