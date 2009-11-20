package calhoun.analysis.crf;

import java.io.Serializable;
import java.util.List;

import calhoun.analysis.crf.io.TrainingSequence;

/** evaluates CRF feature functions on the model.  This is the base class for all feature managers defined within the CRF.<p>
 * <h2>Feature managers versus features</h2>
 * Each feature manager can manage zero, one, or more individual features in the CRF.  This distinction between feature managers and
 * individual features allows a single <code>FeatureManager</code> object to provide evaluations for many related features.  For natural language
 * applications, a single feature manager might hold thousands of features, each corresponding to the presence of a single word in its dictionary.  
 * In the special case of a feature manager that contains only contraints and never has non-zero evaluations, the feature manager may contain zero features.
 * <h2>Feature manager lifecycle</h2>
 * There are three basic stages to the <code>FeatureManager</code> lifecycle.
 * <ul>
 * <li><b>Configuration</b> - Instantiating the feature manager from XML, setting the input component and other options.
 * <li><b>Training</b> - Parameterizing the features on a set of training data.
 * <li><b>Evaluation</b> - Evaluating the features during training or inference.
 * </ul>
 * <h3>Configuration</h3>
 * <code>FeatureManager</code>s are normally specified as part of a model in an XML model file and are instantiated and configured 
 * by the engine when the model file is read.  It is possible to add arbitrary configuration properties to any <code>FeatureManager</code>
 * by adding public get and set methods for the properties.  For more detail on this, see <a href="http://www.springframework.org/documentation">the Spring documentation</a>.
 * <p>
 * The one piece of configuration common to all <code>FeatureManager</code>s is the <code>inputComponent</code> property.  This property is
 * used when the input data contains multiple different sources of input, and the <code>FeatureManager</code> computes its evaluations off of 
 * only one of them.  In this case, each input source is given a different component name.  The <code>inputComponent</code> property of the
 * <code>FeatureManager</code> is then set to the name of the input source it uses for evaluation (usually this is done within the XML configuration).
 * Once this configuration is complete, the Conrad engine will ensure that the feature manager is handed only the correct component during evaluation.
 * This enables modular feature managers.  A feature manager can be written to
 * assume a very specific set of input data.  It can then be combined together in a model with feature managers that expect different input data and
 * provided the input components are assigned correctly each feature will behave as expected without any modifications.       
 * <h3>Training</h3>
 * <h4>Setting evaluation parameters during training</h4>
 * The training stage of the lifecycle allows the feature manager to examine a full set of training data in case any of its evaluations require parameters
 * that need to be learned.  Although the feature manager is required to implement the <code>train</code> function, it may not parameterize its evalautions.
 * Note that training the features is different from the training the feature weights, which inolves an optimization involving all features.  As an example 
 * of a feature manager which uses training, the DenseWeightedEdge feature examines the training data to compute transition log probabilities for each edge in
 * the model.  Then during the evaluation phase, this log probability is the value returned for each edge.
 * <h4>Determining the number of features</h4>
 * During training, a feature manager must fix the number of features it's managing.  For some feature managers this number may be fixed beforehand,
 * but for others it may be not be determined until the training data is examined.  An example is a feature manager that creates a dictionary on text and has one feature for each word in the input
 * data.  The exact number of features might not be determined until the training data is examined and the number of unique words is determined.  After the train
 * function has been called, it is legal for the engine to call {@link #getNumFeatures()}.
 * <p>
 * Additionally, during the call to <code>train</code> the engine must assign feature indices to each feature.  Since many feature managers may be combined together
 * and each feature requires a unique index, each feature manager is given a range of number to use for indices.  The range begins with the <code>startIx</codeIx>
 * parameter that is passed into train and includes as many consecutive integers as the feature manager has features.  Therefore, a feature manager with 3 features
 * that receives a <code>startIx</code> of 5 will assign feature ids 5, 6, and 7 to its features.  The next feature manager will then be given 8 as a <code>startIx</code>
 * <h3>Evaluation</h3>
 * Once a feature manager is configured and trained it is available to the Conrad engine for evaluation.  Each <code>FeatureManager</code> must implement one or
 * more of the <code>FeatureEvaluation</code> interfaces.  Each of these interfaces defines an evaluate function.  The most general evaluate function is a real-valued 
 * function which is dependent on the entire input sequence, a start and end position in the hidden sequence, a hidden state, and a previous hidden state.  
 * During training and inference, each feature must be evaluated for each state at every position.  Therefore <code>evaluate</code> will be called many times
 * during each of these processes.
 * <p>
 * When evaluate is called, the <code>FeatureManager</code> is responsible for computing a value for each feature it manages.  For each non-zero evaluation, the
 * feature manager must call {@link calhoun.analysis.crf.FeatureList#addFeature(int, double)} with the feature index and the evaluated value.  Additionally, feature
 * managers can also declare a particular state, transition, or segment invalid.  In that case, the engine will remove all paths containing that state, transition or segment
 * from consideration.
 * <h2>Other notes</h2> 
 * <h3>Serialization</h3>
 * All features should be serializable.  Conrad uses Java serialization to save its models, and so a feature must be serializable if it is to be 
 * saved as part of a trained model.  The advantage of using serialization is that all parameterization of models will be automatically stored with the model,
 * making it easy to build and incorporate new features into a model.  As a side effect of this, features should not store unnecessarily large amounts of data, such
 * as a full copy of the training data, as this will be stored in the serialized model and read in again for inference.  This will slow the reading of the model and
 * add to memory overhead during inference.
 * <h3>Performance</h3>
 * Conrad does extensive caching and attempts to make the minimum number of calls possible to the feature evaluate functions.  Thus, during training iterations, the engine
 * is generally working with cached versions of the feature evaluations and the evaluate functions are not called.  Therefore, extensive optimization of <code>FeatureManager</code>
 * code is not usually necessary to achieve adequate training performance.  For inference however, the bottleneck is usually in the evaluation of the features, since 
 * inference does not require iterative evaluation of the features. 
 * <h3>Choosing a derived class</h3> 
 * The are several <code>FeatureEvaluation</code> interfaces to choose from.  The various interfaces each take different subsets of the most general feature function
 * signature.  It is important to choose the correct interface when implementating a feature manager.  Choosing a feature manager with extra parameters can lead to 
 * inefficient caching and slow performance.  Choosing one with too few parameters can lead to incorrect results.
 * @param <InputType> the InputSequence type this FeatureManager expects.
 */
public interface FeatureManager<InputType> extends Serializable {

	/** caching strategy that the {@link calhoun.analysis.crf.solver.CacheProcessor} should use to cache values for this feature.  
	 * This is only a hint, the cache processor is not required to use this (or any) caching strategy.  This base class defaults
	 * to the UNSPECIFIED cache strategy.
	 * 
	 * @return cache strategy specification appropriate for this feature.
	 */ 
	CacheStrategySpec getCacheStrategy();
	
	/** Returns a human identifiable name for the feature referenced by a given index.  Used for display purposes only.
	 * @param featureIndex	the index of this feature
	 * @return the human readable name of this feature
	 * */
	String getFeatureName(int featureIndex);

	/** Specifies the particular component of the input which this feature manager uses, if the input data is a composite input which has multiple components.  
	 * May be <code>null</code> if the input is not a composite or if the feature manager has access to all
	 * of the input.  See the <i>Conrad User's Guide</i> for <a href="golem:8080/display/Conrad/Specifying_Input_Data">how to set up input data</a>.
	 * @return string name of the input component this feature should look at, or null if the feature has access to all inputs. */  
	String getInputComponent();

	/** Sets which particular component of a CompositeInput this FeatureManager has access to.
	 *	@param name name of the component of the CompositeInput this FeatureManager has access to.  Usually set automatically during configuration. 
	 */
	void setInputComponent(String name);

	/** Provides access to the entire training set so that FeatureManager can compute global properties and assign feature indices.  This will be called before
	 * any evaluations are requested.  If the FeatureManager can have a variable number of features, this must be fixed within this method.
	 *  
	 * @param startingIndex 	the feature index of the first feature owned by this FeatureManager.  Each FeatureManager 
	 * must use up consecutive indexes, so the last index used will be startingIndex + numFeatures - 1. 
	 * @param modelInfo 		the model that contains this feature 
	 * @param data 				the full list of training sequences to use to train the feature 
	 * */
	void train(int startingIndex, ModelManager modelInfo, List<? extends TrainingSequence<? extends InputType> > data);

	/** Returns the number of features maintained by this <code>FeatureManager</code>.  This number must be fixed after the call to trainFeatures is complete.
	 * @return number of features managed by this <code>FeatureManager</code>
	 */
	int getNumFeatures();
}
