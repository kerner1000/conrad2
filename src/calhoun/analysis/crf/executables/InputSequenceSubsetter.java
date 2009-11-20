package calhoun.analysis.crf.executables;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.Conrad;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.util.Assert;
import calhoun.util.FileUtil;

/** takes an input sequence and chops it up into smaller pieces.  Useful for generating training and test sets.<p>
The proper usage is [configfile] [inputfile] [regionsfile] [pad] [outputfile] [flagForceGenic]
 */
public class InputSequenceSubsetter {
	private static final Log log = LogFactory.getLog(InputSequenceSubsetter.class);

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws Exception {
		if (args.length != 6) {
			System.out.println("The proper usage is [configfile] [inputfile] [regionsfile] [pad] [outputfile] [flagForceGenic]");
			System.out.println("  Note: if you don't care whether the region(s) are gene(s) or not, then set pad=0, flag=0.");
			System.out.println("  A regions file is tab delimited with [seqname] [start] [end], one line per region");
			Assert.a(false);
		}
		String configFile  = args[0];  // "test/input/test_subsetting_configfile.txt";
		String inputFile   = args[1];  // "test/input/test_subsetting_inputfile.txt";  
		String regionsFile = args[2];  // "test/input/test_subsetting_regions.txt"
		int pad            = Integer.parseInt(args[3]);
		String outputFile  = args[4];  // "test/output/test_subsetting_output.txt";
		int flagForceGenic = Integer.parseInt(args[5]);
			
		
		// Define a model manager from file
		Conrad c = new Conrad(configFile);
		ModelManager cm= c.getModel();
		
		// Read a list of inputsequences from file, of the type expected by model manager cm
		Iterator<? extends TrainingSequence<?>> iter = c.getInputHandler().readTrainingData(inputFile).iterator();
		
		String[][] regions  = FileUtil.readFlatFile(regionsFile);
		System.out.println("Number of regions in regions file is " + regions.length);
		
		List<TrainingSequence<?>> s = new ArrayList();
		
		while(iter.hasNext()) {
			TrainingSequence<?> t = iter.next();
			
			String targetseqname = (String) t.getInputSequence().getComponent("name").getX(0);
			
			for (int j=0; j<regions.length; j++) {
				Assert.a( regions[j].length == 3);		
				
				String seqname = regions[j][0];
				
				if (targetseqname.equals(seqname)) {
					int start      =  Integer.parseInt(regions[j][1])-pad;
					int end        =  Integer.parseInt(regions[j][2])+pad;
					log.debug("Region: "+seqname+": "+start+"-"+end);

					if ( (start < 1) || (start>end) || (end>t.length()) ) { 
						log.debug("Skipping Region: "+seqname+": "+start+"-"+end);
						continue;
					}
					
					TrainingSequence u = t.subSequence(start,end); 
					int	len = u.length();
					System.out.println("OK so far and length = " + len);
					if (flagForceGenic > 0) {
						for (int k=0; k<pad; k++) {
							if (u.getY(k) != cm.getStateIndex("intergenic")) { continue; }
						}
						if (u.getY(pad+3) == cm.getStateIndex("intergenic")) { continue; }
						if (u.getY(len-pad-4) == cm.getStateIndex("intergenic")) { continue; }
						for (int k=len-pad; k<len; k++) {
							if (u.getY(k) != cm.getStateIndex("intergenic")) { continue; };						
						}
					}
					s.add(u);
				}
			}
		}
		
		// Write these subsetted InputSequences to file; they should be of exactly the same type as before.
		System.out.println("Will write " + s.size() + " training sequences to file");
		c.getInputHandler().writeTrainingData(outputFile, s);
	}
}
