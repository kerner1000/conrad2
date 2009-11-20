package calhoun.analysis.crf;

import java.io.IOException;
import java.io.Serializable;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.springframework.beans.factory.xml.XmlBeanDefinitionReader;
import org.springframework.context.ApplicationContext;
import org.springframework.context.support.GenericApplicationContext;
import org.springframework.core.io.ByteArrayResource;

import calhoun.analysis.crf.CRFInference.InferenceResult;
import calhoun.analysis.crf.io.InputHandler;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.analysis.crf.io.OutputHandler;
import calhoun.analysis.crf.io.OutputHandlerGeneCallPredict;
import calhoun.analysis.crf.io.OutputHandlerGeneCallStats;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.analysis.crf.io.CompositeInput.LegacyInputHandler;
import calhoun.util.ErrorException;
import calhoun.util.FileUtil;

/** the central class for the Conrad engine.  Has a <code>main</code> function for calling Conrad
 * from the command line and a high-level public interface for programmatic operations.  This class
 * is mostly just a container which delegates the real work to the various objects set up in the
 * configuration.
 * */
public class Conrad implements Serializable {
	private static final long serialVersionUID = -5964610632818921236L;
	private static final Log log = LogFactory.getLog(Conrad.class);
	
	byte[] configXml;
	ModelManager model;
	InputHandler inputHandler;
	OutputHandler outputHandler;
	transient CRFInference inference;
	transient CRFTraining optimizer;
	
	double[] weights = null;
	double trainingTime = 0;
	
	/** Command line entry point for running CRFs.  Can train, test, or predict depending on arguments.  
	 * @param args list of command line arguments.  See usage for details. */
	public static void main(String[] args) throws Exception {
		if(args.length != 4) {
			usage();
		}
		if(args[0].startsWith("train")) {
			Conrad crf;
			if(args[0].equals("trainWeights")) {
				crf = Conrad.read(args[1]);
				crf.trainWeights(args[2]);
			}
			else {
				crf = new Conrad(args[1]);
				if(args[0].equals("trainFeatures")) {
					crf.trainFeatures(args[2]);
				}
				else {
					crf.train(args[2]);			// train features and weights
				}
			}
			crf.write(args[3]);
		}
		else if(args[0].equals("test")) {
			Conrad crf = Conrad.read(args[1]);
			crf.initSolver();
			crf.test(args[2], args[3]);
		}
		else if(args[0].equals("predict")) {
			Conrad crf = Conrad.read(args[1]);
			OutputHandlerGeneCallPredict predictOutputHandler = new OutputHandlerGeneCallPredict();
			predictOutputHandler.setWriteTrainingData(false);
			predictOutputHandler.setManager(crf.getModel());
			predictOutputHandler.setInputHandler(crf.getInputHandler());
			crf.setOutputHandler(predictOutputHandler);
			crf.initSolver();
			crf.testWithoutAnswers(args[2], args[3]);
		}
		else {
			usage();
		}
	}

	/** creates a Conrad engine with no configuration.  All configuration must be done programmatically. 
	 */
	public Conrad() {
	}
	
	/** creates a Conrad engine based on configuration information from an XML model file.
	 * @param configFile string filename of the XML model file  */
	@SuppressWarnings("deprecation")
	public Conrad(String configFile) {
		try {
			configXml = FileUtil.readFileAsBytes(configFile);
		}
		catch(IOException ex) {
			throw new RuntimeException("Failed to load config file: "+configFile, ex);
		}
		ApplicationContext ctx = initSolver();
		model = (ModelManager) ctx.getBean("model");
		
		// For backwards compatibility, if we don't have an inputHandler but we do have an input format, use a legacy handler.
		if(ctx.containsBean("inputFormat") && !ctx.containsBean("inputHandler")) {
			inputHandler = new LegacyInputHandler(ctx.getBean("inputFormat"));
			
			if(ctx.containsBean("outputHandler")) {
				outputHandler = (OutputHandler)ctx.getBean("outputHandler");				
			} else {
				OutputHandlerGeneCallStats legacyOutputHandler = new OutputHandlerGeneCallStats();
				legacyOutputHandler.setWriteTrainingData(true);
				legacyOutputHandler.setManager(model);
				legacyOutputHandler.setInputHandler(inputHandler);	
				outputHandler = legacyOutputHandler;
			}			
		}
		else {
			inputHandler = (InputHandler) ctx.getBean("inputHandler");
			outputHandler = (OutputHandler) ctx.getBean("outputHandler");
		}
	}

	/**  writes this Conrad engine to a file.  This is most often used to save a trained model file.
	 * 
	 * @param filename string name of the file that will contain the serialized model.
	 * @throws IOException if a problem occurs writing to the file
	 */
	public void write(String filename) throws IOException {
		FileUtil.writeObject(filename, this);
	}

	/** read in a Conrad engine from a file.  This file must have previously been created by calling {@link #write}
	 * 
	 * @param filename string name of the file containing the model.
	 * @return the Conrad engine which has been read from the file
	 * @throws IOException if there is a problem reading the file
	 */
	public static Conrad read(String filename) throws IOException {
		try {
			Conrad ret = (Conrad) FileUtil.readObject(filename);
			ret.initSolver();
			return ret;
		} catch (ClassNotFoundException ex) {
			throw new ErrorException(ex);
		}
	}

	/** fully trains this Conrad engine with this training data.  The training data is specified as a string location,
	 * which the configured {@link InputHandler } is responsible for converting into a list of training sequences
	 * of the appropriate type.
	 * 
	 * @param location string location of the data.  The exact meaning will be determined by the InputHandler.
	 * @throws IOException if there is a problem reading the training data.
	 */
	public void train(String location) throws IOException {
		List<? extends TrainingSequence<?>> t = inputHandler.readTrainingData(location, false);
		train(t);
	}

	
	/** fully trains this Conrad engine with this training data.  The training data is specified as a list of 
	 * training sequences, and no DataInputHandler is used.
	 * @param data a list of training sequences to use for training
	 */
	public void train(List<? extends TrainingSequence<?>> data) {
		trainFeatures(data);
		trainWeights(data);
	}
	
	/** trains only the features in the current model with this training data.  {@link FeatureManager#train } is called
	 * for each feature in the model, but no optimization is performed and no feature weights are set.  This allows the
	 * features themselves to be parameterized on one set of training data, while using a different set for optimizing the
	 * feature weights.<p>
	 * The training data is specified as a string location,
	 * which the configured {@link InputHandler } is responsible for converting into a list of training sequences
	 * of the appropriate type.
	 * 
	 * @param location string location of the data.  The exact meaning will be determined by the InputHandler.
	 * @throws IOException if there is a problem reading the training data
	 */
	public void trainFeatures(String location) throws IOException {
		List<? extends TrainingSequence<?>> t = inputHandler.readTrainingData(location, false);
		trainFeatures(t);
	}

	/** trains only the features in the current model with this training data.  {@link FeatureManager#train } is called
	 * for each feature in the model, but no optimization is performed and no feature weights are set.  This allows the
	 * features themselves to be parameterized on one set of training data, while using a different set for optimizing the
	 * feature weights.<p>
	 * The training data is specified as a list of 
	 * training sequences, and no DataInputHandler is used.	 
	 * @param data a list of training sequences to use for training
	 */
	public void trainFeatures(List<? extends TrainingSequence<?>> data) {
		print("Training features");
		double timer = System.currentTimeMillis();

		// Train features
		model.train(0, model, data);
		if (log.isDebugEnabled()) {
			log.debug("Features:");
			for (int i = 0; i < model.getNumFeatures(); ++i) {
				log.debug(model.getFeatureName(i));
			}
		}
		
		trainingTime = (System.currentTimeMillis() - timer)/1000;
		print("Trained in "+trainingTime+" seconds.");
	}

	/** optimizes the feature weights for the current model with this training data.  Assumes that {@link #trainFeatures } 
	 * has already been called to train the individual features.  <p>
	 * The training data is specified as a list of  training sequences, and no DataInputHandler is used.	  
	 * @param location string location of the data.  The exact meaning will be determined by the InputHandler.
	 * @throws IOException if there is a problem reading the training data
	 */
	public void trainWeights(String location) throws IOException {
		List<? extends TrainingSequence<?>> t = inputHandler.readTrainingData(location, false);
		trainWeights(t);
	}

	/** optimizes the feature weights for the current model with this training data.  Assumes that {@link #trainFeatures } 
	 * has already been called to train the individual features.  <p>
	 * The training data is specified as a list of  training sequences, and no DataInputHandler is used.	 
	 * @param data a list of training sequences to use for training
	 */
	public void trainWeights(List<? extends TrainingSequence<?>> data) {
		print("Training weights");
		double timer = System.currentTimeMillis();

		// Train weights
		weights = optimizer.optimize(model, data);

		timer =  (System.currentTimeMillis() - timer)/1000;
		trainingTime += timer;
		print("Trained weights in "+timer+" seconds.  "+trainingTime+" total.");
	}

	/** runs a trained model against a set of input data with known results and evaluates the performance.  Assumes that
	 * {@link #train} has already been called to train the model.  For convenience, the data is passed in as a training set,
	 * although the model is not trained.  The input is used to create a set of predictions and then those predictions are
	 * compared against the expected outputs.  The result of the prediction is passed to the output handler which can compare
	 * the predicted versus the expected values
	 * @param inputLocation string location of the data.  The exact meaning will be determined by the InputHandler.
	 * @throws IOException if there is a problem reading the training data
	 */
	public void test(String inputLocation) throws IOException {
		List<? extends TrainingSequence<?>> t = inputHandler.readTrainingData(inputLocation, false);
		test(t);
	}

	public void test(String inputLocation, String outputLocation) throws IOException {
		List<? extends TrainingSequence<?>> t = inputHandler.readTrainingData(inputLocation, false);
		test(t, outputLocation);
	}

	public void test(List<? extends TrainingSequence<?>> data) throws IOException {
		test(data, null);
	}
	
	public void testWithoutAnswers(String inputLocation, String outputLocation) throws IOException {
		List<? extends TrainingSequence<?>> t = inputHandler.readTrainingData(inputLocation, true);
		test(t, outputLocation);
	}
	
	/** runs a trained model against a set of input data with known results and evaluates the performance.  Assumes that
	 * {@link #train} has already been called to train the model.  For convenience, the data is passed in as a training set,
	 * although the model is not trained.  The input is used to create a set of predictions and then those predictions are
	 * compared against the expected outputs.  The result of the prediction is passed to the output handler which can compare
	 * the predicted versus the expected values
	 * @param data a list of training sequences to use for training
	 */
	public void test(List<? extends TrainingSequence<?>> data, String location) throws IOException {
		print("Beginning test");
		printWeights();
		outputHandler.setOutputLocation(location);
		for (TrainingSequence dr : data) {
			InferenceResult predictedHiddenSequence = predict(dr);
			outputHandler.writeTestOutput(dr.getInputSequence(), dr.getY(), predictedHiddenSequence.hiddenStates);
		}
		print("Testing complete");
		outputHandler.outputComplete();
	}

	/** preforms inference on the input sequence and determines the best labeling for the sequence using
	 * the configured inference algorithm.
	 * @param data the input sequence the engine will use for inference
	 * @return an inference result containing the predicted hidden states
	 */
	public InferenceResult predict(InputSequence data) {
		return inference.predict(model, data, weights);
	}

	/** sets feature weights.  Usually these weights are determined during the training process, but they can be set directly.
	 * @param weights an array of doubles containing one weight for each feature.
	 */
	public void setWeights(double[] weights) {
		this.weights = weights;
	}
	
	/** looks up a feature's name given it's index
	 * @param index index of the feature
	 * @return name of the feature
	 */
	public String getFeatureName(int index) {
		return model.getFeatureName(index);
	}

	/** returns the the total number of seconds used in training.  This is the sum of the time to train the features and the time to train the weights.
	 * This is set when each phase of the training (features & weights) is completed.
	 * @return total training time
	 */
	public double getTrainingTime() {
		return trainingTime;
	}

	/** returns the number of individual features in the model.  This may differ from the number of <code>FeatureManager</code>s because each
	 * <code>FeatureManager</code> may have 0, 1, or many features associated with it.
	 * @return total number of features in the model.
	 */
	public int getNumFeatures() {
		return model.getNumFeatures();
	}

	/** returns the number of hidden states in the model
	 * @return number of hidden states in the model
	 */
	public int getNumStates() {
		return model.getNumStates();
	}

	/** looks up the name of a state given it's index.
	 * @return string name of the state with this index.
	 */
	public String getStateName(int state) {
		return model.getStateName(state);
	}

	/** returns the configured ModelManager object.
	 * @return the model manager which contains the features and hidden state configuration
	 */
	public ModelManager getModel() {
		return model;
	}
	
	/** returns the configured numerical optimizer which will be used to
	 * select the optimal feature weights during training.
	 * @return the configured objective function gradient
	 */
	public CRFTraining getOptimizer() {
		return optimizer;
	}
	
	/** returns the feature weights.  These will be valid once the modle is trained.
	 * @return an array of doubles containing the weight for each feature.  It will be the
	 * same length as returned by {@link #getNumFeatures()}
	 */
	public double[] getWeights() {
		return weights;
	}

	/** returns the configured inference algorithm which will be used 
	 * to predict hidden states for new inputs once the model is trained.
	 * @return the configured inference algorithm
	 */
	public CRFInference getInference() {
		return inference;
	}

	/** sets the inference algorithm.  Called automatically during configuration. */
	public void setInference(CRFInference inference) {
		this.inference = inference;
	}

	/** sets the model.  Called automatically during configuration. */
	public void setModel(ModelManager model) {
		this.model = model;
	}

	/** sets the numerical optimizer.  Called automatically during configuration. */
	public void setOptimizer(CRFTraining optimizer) {
		this.optimizer = optimizer;
	}

	/** Returns a formatted string listing the weights.  Useful for debugging. 
	 * @return a string containing a human readable list of the feature weights
	 * */
	public String printWeights() {
		StringBuilder sb = new StringBuilder();
		for(int i=0; i<weights.length; ++i) {
			sb.append(String.format("%.5f\t%s\n", weights[i], getFeatureName(i)));
		}
		return sb.toString();
	}

	private void print(String msg) {
		System.out.println(msg);
	}

	private ApplicationContext initSolver() {
		GenericApplicationContext ctx = new GenericApplicationContext();
		XmlBeanDefinitionReader xmlReader = new XmlBeanDefinitionReader(ctx);
		xmlReader.loadBeanDefinitions(new ByteArrayResource(configXml));
		ctx.refresh();
		inference = (CRFInference) ctx.getBean("inference");
		optimizer = (CRFTraining) ctx.getBean("optimizer");
		return ctx;
	}
	
	private static void usage() {
		System.out.println("       Conrad train(Features) configFile data modelFile");
		System.out.println(" or    Conrad trainWeights modelFileIn data modelFileOut");
		System.out.println(" or    Conrad test modelFile inputData outputData");
		System.out.println(" or    Conrad predict modelFile inputData outputData");
		System.exit(-1);
	}

	/** returns the configured input handler.
	 * @return the input handler for this model
	 */
	public InputHandler getInputHandler() {
		return inputHandler;
	}

	/** sets the configured input handler.  Must be set before any train or test methods are called.  Usually called during config based on setup in the XML file.
	 * @param inputHandler the input handler for this model
	 */
	public void setInputHandler(InputHandler inputHandler) {
		this.inputHandler = inputHandler;
	}

	/** gets the configured output handler.  Must be set before any test methods are called.
	 * @return Returns the outputHandler.
	 */
	public OutputHandler getOutputHandler() {
		return outputHandler;
	}

	/** sets the configured output handler.  WillMust be set before any test methods are called.
	 * @param outputHandler the output handler for this model
	 */
	public void setOutputHandler(OutputHandler outputHandler) {
		this.outputHandler = outputHandler;
	}
}
