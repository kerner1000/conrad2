package calhoun.analysis.crf.io;

import java.io.IOException;

/** writes output to a file using a {@link TrainingSequenceIO}.  The name of the output file will be computed from the output
 * location by using the configured {@link FilenameMapper}.  If no mapper is configured, the output location will be used as
 * the filename.
 */
public class OutputHandlerWriter implements OutputHandler {
	private static final long serialVersionUID = -3039234501926777322L;

	FilenameMapper filenameMapper;
	TrainingSequenceIO resultHandler;
	String location;
	
	public void setOutputLocation(String location) {
		this.location = location;
	}
	
	public void writeOutput(InputSequence<?> sequence, int[] hiddenStates) throws IOException {
		throw new UnsupportedOperationException();
	}

	public void writeTestOutput(InputSequence<?> sequence, int[] truePath, int[] hiddenStates) throws IOException {
	}

	public void outputComplete() throws IOException {
	}

	/** gets the mapper which will be used to generate the output filename based on the output location.
	 * @return the configured filenameMapper
	 */
	public FilenameMapper getFilenameMapper() {
		return filenameMapper;
	}

	/** sets the mapper which will be used to generate the output filename based on the output location.
	 * @param filenameMapper the filenameMapper used to generate the output filename based on the output location.
	 */
	public void setFilenameMapper(FilenameMapper filenameMapper) {
		this.filenameMapper = filenameMapper;
	}

	/** gets the result handler which will be used to handle the output sequence
	 * @return the result handler which will be used to handle the output sequence
	 */
	public TrainingSequenceIO getResultHandler() {
		return resultHandler;
	}

	/** sets the result handler which will be used to handle the output sequence
	 * @param resultHandler the result handler which will be used to handle the output sequence
	 */
	public void setResultHandler(TrainingSequenceIO resultHandler) {
		this.resultHandler = resultHandler;
	}
}
