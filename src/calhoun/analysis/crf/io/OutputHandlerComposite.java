package calhoun.analysis.crf.io;

import java.io.IOException;
import java.util.List;

/** an output handler that allows multiple output handlers to be configured.  Each result is handed to each 
 * output handler in the list.
 */
public class OutputHandlerComposite implements OutputHandler {
	private static final long serialVersionUID = 7923868199821973396L;

	List<OutputHandler> handlers;
	
	public void setOutputLocation(String location) {
		for(OutputHandler handler : handlers) {
			handler.setOutputLocation(location);
		}
	}
	
	public void writeOutput(InputSequence<?> sequence, int[] hiddenStates) throws IOException {
		for(OutputHandler handler : handlers) {
			handler.writeOutput(sequence, hiddenStates);
		}
	}

	public void writeTestOutput(InputSequence<?> sequence, int[] truePath, int[] hiddenStates) throws IOException {
		for(OutputHandler handler : handlers) {
			handler.writeTestOutput(sequence, truePath, hiddenStates);
		}
	}

	public void outputComplete() throws IOException {
		for(OutputHandler handler : handlers) {
			handler.outputComplete();
		}
	}

	/** retursn the configured list of output handler.
	 * @return returns the output handlers
	 */
	public List<OutputHandler> getHandlers() {
		return handlers;
	}

	/** sets the list of output handlers 
	 * @param handlers a list of output handlers.  The output handlers will be called in order for each prediction or test result.
	 */
	public void setHandlers(List<OutputHandler> handlers) {
		this.handlers = handlers;
	}
}
