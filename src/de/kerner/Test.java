package de.kerner;
import java.io.IOException;
import java.util.List;

import calhoun.analysis.crf.Conrad;
import calhoun.analysis.crf.io.InputHandler;
import calhoun.analysis.crf.io.OutputHandler;
import calhoun.analysis.crf.io.TrainingSequence;


public class Test {
	
	private final static String TRAIN = "train";
	private final static String MODEL = "models/singleSpecies.xml";
	private final static String PATH_IN = "/home/proj/kerner/diplom/f.graminearum/";
	private final static String FILE_OUT = "/home/proj/kerner/diplom/f.graminearum/fg_3_supercontigs_train.bin";
	
//	private final InputHandler inputHandler = new InputHandlerImpl();
	//private final OutputHandler outputHandler;
	
	public static void main(String args[]){
		final Test w = new Test();
		w.runTraining();
	}
	
	private void runTraining(){
		final Conrad crf = new Conrad(MODEL);
//		crf.setInputHandler(inputHandler);
		//crf.setOutputHandler(outputHandler);
		try {
			crf.trainFeatures(PATH_IN);
			crf.trainWeights(PATH_IN);
			crf.write(FILE_OUT);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
	}

}
