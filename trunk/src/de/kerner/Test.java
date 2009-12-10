package de.kerner;
import java.io.IOException;
import java.util.List;

import calhoun.analysis.crf.Conrad;
import calhoun.analysis.crf.io.InputHandler;
import calhoun.analysis.crf.io.OutputHandler;
import calhoun.analysis.crf.io.TrainingSequence;


public class Test {
	
	public static void main(String args[]){
		predict();
	}
	
	static void predict(){
		final String s1 = "predict";
		final String s2 = "/home/proj/kerner/diplom/anna/annaPredictThaliana/data/conrad/trainingFile.bin";
		final String s3 = "/home/proj/kerner/diplom/anna/annaPredictThaliana/data/conrad";
		final String s4 = "predicted";
		
		try {
			Conrad.main(new String[]{s1, s2, s3, s4});
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
//	private void runTraining(){
//		final Conrad crf = new Conrad(MODEL);
////		crf.setInputHandler(inputHandler);
//		//crf.setOutputHandler(outputHandler);
//		try {
//			crf.trainFeatures(PATH_IN);
//			crf.trainWeights(PATH_IN);
//			crf.write(FILE_OUT);
//		} catch (Exception e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
		
		
//	}

}
