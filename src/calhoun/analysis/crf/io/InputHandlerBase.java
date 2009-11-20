package calhoun.analysis.crf.io;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

/** Utility class to abstract out some of the common behavior for {@link InputHandler}s.
 */
public abstract class InputHandlerBase implements InputHandler {
	private static final long serialVersionUID = 4002620265503458109L;

	Iterator<? extends InputSequence<?>> createCompositeInput(List<Map<String, InputSequence<?>>> inputs) {
		// Create a composite input
		List<InputSequenceComposite> ret = new ArrayList();
		for(Map<String, InputSequence<?>> input : inputs) {
			ret.add(new InputSequenceComposite(input));
		}
		return ret.iterator();
		
	}
	
	List<? extends TrainingSequence<?>> readTrainingData(String location, String trainingLocation, TrainingSequenceIO hiddenStateReader, boolean predict) throws IOException {
		Iterator<? extends InputSequence<Map<String, Object>>> inputIter = (Iterator<? extends InputSequence<Map<String, Object>>>) readInputData(location);
		
		List<TrainingSequence<Map<String, Object>>> ret = new ArrayList<TrainingSequence<Map<String, Object>>>();
		while(inputIter.hasNext()) {
			ret.add(new TrainingSequence(inputIter.next()));
		}
		
		TrainingSequenceIO readerToUse;
		if (predict) {
			readerToUse = new AllIntergenicHiddenStateReader();
		} else {
			readerToUse = hiddenStateReader;
		}
		
		readerToUse.readTrainingSequences(trainingLocation, ret);
		return ret;
	}

}
