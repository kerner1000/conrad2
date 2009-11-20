package calhoun.analysis.crf.solver;

import calhoun.analysis.crf.LocalPathSimilarityScore;
import calhoun.analysis.crf.solver.semimarkov.CleanLocalScoreSemiMarkovGradient;

/** computes an objective function which is the expected value of a local path similarity score on a 
 * semi-Markov model.  Requires a {@link CacheProcessor} and a {@link LocalPathSimilarityScore} to be configured.<p>
 * <h2>Debugging output</h2>
 * To get a better understanding of what the objective function is doing, several different properties can be set that
 * cause the objective function to write out trace files showing its calculations during training.  Usually when turning
 * these options on, you should set <code>maxIters = 1</code> and <code>requireConvergence = false</code> in your optimizer
 * to do only a single training iteration, possibly setting the starts to some predetermined value.  Each of these
 * properties can be configured with a filename and each time {@link #apply} is called, the file will be overwritten with 
 * data from the current call.  The logging options are:
 *
 * <ul>
 * <li> <b><code>alphaFile</code></b> - computation of alpha values for Markov states, includes all nodes and edges.
 * <li> <b><code>alphaLengthFile</code></b> - computation of alpha values for semi-Markov states , includes all segments
 * <li> <b><code>betaLengthFile</code></b> - computation of beta values for semi-Markov states , includes all segments
 * <li> <b><code>expectFile</code></b> - computation of expected values for each Markov feature 
 * <li> <b><code>expectLengthFile</code></b> - computation of expected values for each semi-Markov feature  
 * <li> <b><code>nodeMarginalFile</code></b> - computation of marginal probability of each state at each position 
 * </ul>

 * <h4>Implementation Notes</h4>
 * The general normalization scheme works as follows. When updating alpha values in the forward pass we compute segments
 * of length 1 first and then work backwards.
 * <p>
 * Instead of always normalizing to 1 we discretize the normalization. We choose an arbitrary normalization factor w,
 * such as 50. The normalization factor at any position is then an integer v, and all entries at that position are
 * alpha[y]*e^(v*w).
 * <p>
 * The normalization can be computed at any position from 1) Elements of the alpha array are summed s 2) v = log(s)/w.
 * By integer division v will always be an appropriate normalizer. It may be positive or negative. 3) All elements of
 * the array are divided by e^(v*w)
 * 
 */
public class MaximumExpectedAccuracySemiMarkovGradient extends CleanLocalScoreSemiMarkovGradient {
}