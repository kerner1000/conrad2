package calhoun.analysis.crf.solver.semimarkov;

import java.util.Arrays;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.solver.LogFiles;
import calhoun.analysis.crf.solver.LookbackBuffer;
import calhoun.analysis.crf.solver.CacheProcessor.FeatureEvaluation;
import calhoun.analysis.crf.solver.CacheProcessor.LengthFeatureEvaluation;
import calhoun.analysis.crf.solver.CacheProcessor.SolverSetup;
import calhoun.analysis.crf.solver.CacheProcessor.StatePotentials;
import calhoun.util.Assert;
import calhoun.util.ColtUtil;
import calhoun.util.FileUtil;

/** Updates a beta entry with a potential from a length dependent feature */
class BetaLengthFeatureProcessor {
	static final Log log = LogFactory.getLog(CleanMaximumLikelihoodSemiMarkovGradient.class);
	static final boolean debug = log.isDebugEnabled();

	final CleanMaximumLikelihoodSemiMarkovGradient parent;
	final SolverSetup modelInfo;
	final LogFiles logs;
	
	int seqOffset;

	double prob;
	
	// This vector tracks the marginal probability of each node, accounting for length dependence
	double[] nodeProb;

	boolean globalArrays = false;
	double[][] betas;
	int[] betaNorms;
	double[][] allNodeProb;
	double[][] allEdgeProb;
	double[][] cumulativeStableLogProb;
	
	public BetaLengthFeatureProcessor(CleanMaximumLikelihoodSemiMarkovGradient parent) {
		this.parent = parent;
		this.modelInfo = parent.modelInfo;
		this.logs = parent.logs;

		nodeProb = new double[modelInfo.nStates];
	}

	public void setGlobalArrays(double[][] betas, int[] betaNorms, double[][] allNodeProb, double[][] allEdgeProb) {
		globalArrays = true;
		this.betas = betas;
		this.betaNorms = betaNorms;
		this.allNodeProb = allNodeProb;
		this.allEdgeProb = allEdgeProb;
	}
	
	/**
	 * Computing the betas is a bit tricky. The problem is that our cache associates each explicit
	 * length feature with the position at which it ends. For the betas, we are going backwards, and so we need to
	 * look at the positions where the length based feature starts.
	 * 
	 * The way we handle this is by evaluating the length based features as they arise and incrementally adding
	 * their potential values into the beta positions as we go along. Therefore we aren't always filling one beta
	 * matrix at a time, we fill any beta entry that starts a state which ends the current position as we move
	 * backwards through the sequence.
	 * 
	 * In order to do this we need to keep mi matrices and arrays of the potentials of staying in a stable state for
	 * the entire lookback period. These matrices let us compute betas anywhere we want.
	 * 
	 * Unfortunately this process needs more exponentiation than the existing process.
	 * 
	 * Normalization works as follows: First, we need to find the right normalization factor for the value we are
	 * adding. 1) We add all of the exponents we will need. 2) We integer divide by NORM_FACTOR to get the
	 * normalization constant for our beta update 3) The constant*NORM_FACTOR is subtracted from the total exponent
	 * to get the unnormalized exponent 4) The unnormalized exponent is exponentiated and multiplied by the
	 * following beta value 5) We compare the norm constant to the current constant for the beta we are about to
	 * update 5.1) If they are equal we are fine, the values can be added 5.2) If Our new norm constant has larger
	 * abs(), we renormalize the beta vector by dividing all elements by e^(NORM_FACTOR*diff) 5.3) If our new norm
	 * constant has smaller abs(), we renormalize it but dividing by e^(FACTOR*diff) 6) When we finally update the
	 * beta vector at the end we redo the normalization anyway.
	 * 
	 * To optimize this for fastest execution it ends up being one big function. Putting function calls in the inner
	 * loops really slowed things down. I think the optimizer can't do as good a job if there are function calls
	 * there.
	 */
	void computeBetasAndExpectations(int seqNum, int len) {
		// lookback positions from the last base
		// Also initialize the array of previous beta vectors
		seqOffset = modelInfo.seqOffsets[seqNum];
		Arrays.fill(nodeProb, 0.0);

		// Begin by initializing the lookback information for the last position
		// Betas for the last position are all 1.  exp(0)
		LookbackBuffer posLookback = parent.nextBuffer;
		posLookback.clear();
		Arrays.fill(posLookback.beta, 1);
		posLookback.betaNorm = 0;
		Arrays.fill(posLookback.stableState, 0);
		posLookback.pos = len-1;

		// Initialize the stable state vector
		// Need to prime with a stable state of 0 since we don't have the init that we do for alphas.
		double[] lastStableState = posLookback.stableState;

		// Put this position in the lookbackbuffer
		parent.nextBuffer = parent.lookbackBuffer.addFirst(posLookback);

		// Now work through the sequence backwards
		// lastInitPos holds the position of the leftmost position which gets evaluated at the start of the run.
		int lastInitPos = len - 2 - modelInfo.maxLookback;
		int miPos = len-2;
		for (int pos = len - 1; pos >= 0; --pos) {
			// First, update the lookback, which caches mi and stable values if necessary
			while(miPos >= 0 && miPos >= lastInitPos) {
				// Update stable states given the previous transition
				LookbackBuffer newLookback = parent.nextBuffer;
				newLookback.clear();
				newLookback.pos = miPos;
				parent.cacheMi(seqNum, newLookback.mi, lastStableState, newLookback.stableState, miPos+1);
				lastStableState = newLookback.stableState;
				parent.nextBuffer = parent.lookbackBuffer.addFirst(newLookback);
				--miPos;
			}
			--lastInitPos;
			// At this point, miPos contains the leftmost position for which markov features have been evaluated.
			// The lookback buffer at this position has a stableState vector for each length dependent state and
			// an mi matrix computed for that state.
			
			// Check that posLookback is correct
			if(debug) {
				Assert.a(posLookback.pos == pos, "Wrong lookback buffer: was ", posLookback.pos, " should be ", pos);
			}

			// Retrieve the lookback information for the current position.  The lookback buffer runs left to right, so the 
			// current position is not at the beginning.
			LookbackBuffer nextLookback = null;
			if(pos != 0) {
				nextLookback = parent.lookbackBuffer.get(pos - (miPos+1)-1);
				Assert.a(nextLookback.pos == pos-1, "Wrong next lookback buffer: was ", posLookback.pos, " should be ", pos-1);
			}

			// Update beta first for the markov states and then the semi-markov states.  For pos 0 the beta is not updated, but node marginals are still calculated.
			regularBetaUpdate(pos, nextLookback, posLookback.beta, posLookback.betaNorm, posLookback.transitionProb);
			
			// At this point, the nodeProb and edgeProb vectors contain the marginal probabilities of markov nodes and edges into them.
			// The beta is fully updated for the markov states.  Semi markov has not been done.
			lengthBeta(seqNum, pos, miPos, posLookback);

			// For the semi-markov states, calculate self-transition marginals
			// Take the node probabilities, and subtract off each transition
			if(pos > 0) {
				for(StatePotentials lb: modelInfo.statesWithLookback) {
					int state = lb.state;
					int selfTransIndex = modelInfo.selfTransitions[state];
					double selfTransProb = nodeProb[state];
					for(byte pot : lb.potentials) {
						double lbTrans = nextLookback.transitionProb[pot - modelInfo.nStates];
						selfTransProb -= lbTrans;
					}
					Assert.a(nextLookback.transitionProb[selfTransIndex] == 0.0);
					nextLookback.transitionProb[selfTransIndex] = selfTransProb;
				}
			}
			
			verifyMarginals(seqNum, pos, nextLookback);
			
			// At this point the previous beta values are all updated.
			// Update expectations for the transitions we just calculated.
			updateExpectations(seqNum, pos, nextLookback);

			// If we are using global arrays (for local score computation), populate them now
			if(globalArrays) {
				for(int i=0; i<modelInfo.nStates; ++i) {
					allNodeProb[i][pos] = nodeProb[i];
					betas[i][pos] = posLookback.beta[i];
					betaNorms[pos] = posLookback.betaNorm;
				}
				for(int i=0; i<modelInfo.nTransitions; ++i) {
					allEdgeProb[i][pos] = posLookback.transitionProb[i];
				}
			}
			
			// For semi-markov states, intialize the node marginals at the next position
			// This sets nodeProb to the marginal probability of all segments for this state that
			// end after the previous position.  Those that end at the previous position are added in the
			// next call to lengthBeta.
			if(pos > 0) {
				for(StatePotentials lb: modelInfo.statesWithLookback) {
					int selfTransIndex = modelInfo.selfTransitions[lb.state];
					nodeProb[lb.state] = nextLookback.transitionProb[selfTransIndex];
				}
			}
			
			posLookback = nextLookback;
		}
	}

	/**
	 * Does the update of the beta values for all of the regular, non-length dependent states. Also calculates all
	 * of the node and edge probabilities for the non-length nodes and the transitions into them.
	 */
	private void regularBetaUpdate(final int pos, final LookbackBuffer nextLookback, final double[] oldBeta, final int oldNorm, final double[] transitionProb) {
		// If we don't have any markov states, we can skip this.  Importantly, we don't have to worry about the normalization.
		if(modelInfo.statesWithoutLookback.length == 0)
			return;
		
		double[] nodeAlpha = parent.alphas[pos];
		double nodeNorm = CleanMaximumLikelihoodSemiMarkovGradient.exp((parent.alphaNorms[pos] + oldNorm - parent.zNorm) * CleanMaximumLikelihoodSemiMarkovGradient.NORM_FACTOR) * parent.zInv;
		double[] edgeAlpha = null;
		double edgeNorm = Double.NaN;
		double normAdjust = 0.0;
		if(pos > 0) {
			edgeAlpha = parent.alphas[pos-1];

			// Need to deal with possibly different normalization constants.
			if (oldNorm > nextLookback.betaNorm) {
				// We never make the constant smaller, so set the constant for the new factor
				//log.info(String.format("Renormalizing beta[%d] from %d to %d", pos, newNorm, oldNorm));
				CleanMaximumLikelihoodSemiMarkovGradient.renormalize(nextLookback.beta, nextLookback.betaNorm, oldNorm);
				nextLookback.betaNorm = oldNorm;
			} else {
				// The case where beta(pos-1) already has a normalization constant larger than beta(pos)
				// We need to adjust all of the updated beta values.
				normAdjust = (oldNorm - nextLookback.betaNorm) * CleanMaximumLikelihoodSemiMarkovGradient.NORM_FACTOR;
			}

			// We add newNorm here because it ends up cancelling with normAdjust when we calc the edge prob.
			edgeNorm = CleanMaximumLikelihoodSemiMarkovGradient.exp((parent.alphaNorms[pos-1] + nextLookback.betaNorm - parent.zNorm) * CleanMaximumLikelihoodSemiMarkovGradient.NORM_FACTOR) * parent.zInv;
		}

		for(StatePotentials potentials : modelInfo.statesWithoutLookback) {
			byte node = potentials.state;
			double nodePotential = 0.0;
			double betaVal = oldBeta[node];
			if(logs.nodeMarginalWriter != null) {
				FileUtil.safeWrite(logs.nodeMarginalWriter, String.format("NodeMarg[%d][%d] = %f = %f * %f * %f (aN: %d bN: %d zN: %d 1/z: %f)\n", pos, node, nodeAlpha[node] * betaVal * nodeNorm, nodeAlpha[node], betaVal, nodeNorm, this.parent.alphaNorms[pos], oldNorm, this.parent.zNorm, this.parent.zInv));
			}
			nodeProb[node] = nodeAlpha[node] * betaVal * nodeNorm;
			
			// For regular states, we sum edge probabilities to get node probabilities.
			if(pos > 0) {
				for (short potential : potentials.potentials) {
					int trans = potential - modelInfo.nStates;
					double transVal = nextLookback.mi[trans];
					if(!Double.isInfinite(transVal)) {
						// Mi for beta is not exponentiated, so we do it here.
						double potentialValue = CleanMaximumLikelihoodSemiMarkovGradient.exp(nextLookback.mi[trans] + normAdjust);
						nodePotential += potentialValue;

						int from = modelInfo.transitionFrom[trans];
						nextLookback.beta[from] += potentialValue * betaVal;
						nextLookback.transitionProb[trans] = edgeAlpha[from] * potentialValue * betaVal * edgeNorm;
					}
					else {
						nextLookback.transitionProb[trans] = 0.0;
					}
				}
			}
		}

		// Now check to see if this new vector needs normalization again.
		if(pos > 0) {
			try {
				nextLookback.betaNorm += CleanMaximumLikelihoodSemiMarkovGradient.normalize(nextLookback.beta);
			}
			catch(RuntimeException ex) {
				CleanMaximumLikelihoodSemiMarkovGradient.log.warn("Normalization problem at " + pos + " "+ColtUtil.format(nextLookback.beta));
				throw ex;
			}
		}
	}

	/**
	 * Updates an existing beta by adding in: potentialValue - The value of any length-dependent features for this
	 * node f(y, i, d) and edge f(y', y, i, d) stableValue - The value of the non-length dependent node features
	 * summed across the length of this segment mis - The value of the non-length dependent transition from the
	 * previous node to this one f(y', y, i-d)
	 * 
	 * In addition to computing the betas for length dependent features, it updates the probabilities vectors needed
	 * for the feature expectations. These are the marginal at each position or each edge.
	 */
	private void lengthBeta(int seqNum, int pos, int miPos, LookbackBuffer posLookback) {
		double[] beta = posLookback.beta;
		int betaNorm = posLookback.betaNorm;
		double[] lengthStable = posLookback.stableState;
		
		parent.cacheProcessor.evaluateSegmentsEndingAt(seqNum, pos);

		int nSemiMarkovStates = modelInfo.statesWithLookback.length;
		for(int i=0; i<nSemiMarkovStates; ++i) {
			LengthFeatureEvaluation[] lookbacksForState = parent.lengthEvals[i];
			StatePotentials statePotentials = modelInfo.statesWithLookback[i];
			byte toNode = statePotentials.state;

			// If there is no way to end a segment at this end position, give up.
			// Usually this happena t boudnary cases at the edges.
			if(beta[toNode] == 0)
				continue;
			
			int lbArrayIndex=0;
			LengthFeatureEvaluation lengthEval = lookbacksForState[lbArrayIndex];
			int lookback = lengthEval.lookback;
			while(lookback != -1) {
				int prevPos = pos - lookback - 1;
				int lbIndex = prevPos - miPos - 1;

				//log.info("Pos: "+pos+"\t State: "+modelInfo.statesWithLookback[i].state+"\t Lookback: "+lookback);
				LookbackBuffer segBegin = null;
				if(prevPos >= 0) {
					// For speed I hand inline RecyclingBuffer.get
					segBegin = parent.lookbackBuffer.array[(parent.lookbackBuffer.currentStart+lbIndex)%parent.lookbackBuffer.length];
				}
				LookbackBuffer stableBuffer = parent.lookbackBuffer.array[(parent.lookbackBuffer.currentStart+lbIndex+1)%parent.lookbackBuffer.length];

				// Handle evaluation of the node potentials
				double stableValue = stableBuffer.stableState[toNode] - lengthStable[toNode];
				double nodePotential = stableValue;
				
				FeatureEvaluation nodeEvals = lengthEval.nodeEval;
				short[] indices = nodeEvals.index;
				float[] vals = nodeEvals.value;
				int ix = 0;
				short index = indices[ix];
				while(index >= 0) {
					nodePotential += vals[ix] * parent.lambda[index];
					index = indices[++ix];
				}
				if(debug)
					Assert.a(index != Short.MIN_VALUE, "Node lengths should only be returned in the cache if they are valid.  They can be invalid because a node is invalid or a self-transition edge is invalid.");

				if(prevPos < 0) {
					// If this is the first segment, then we don't worry about edges and handle the node directly.
					double expVal = nodePotential + parent.starterAlpha[toNode];
					lengthBetaHandling(seqNum, prevPos, pos, expVal, -1, toNode, 1.0, 0, beta[toNode], betaNorm, nodeEvals);
				}
				else {
					// If this is not the first segment, we need to deal with edges
					FeatureEvaluation[] edgeEvals = lengthEval.edgeEvals;
					int nEdges = statePotentials.potentials.length;
					for(int edgeIx=0; edgeIx < nEdges; ++edgeIx) {
						int potential = statePotentials.potentials[edgeIx];
						int trans = potential - modelInfo.nStates;
						int fromNode = modelInfo.transitionFrom[trans];
						// Skip semi-Markov self transitions
						if(fromNode == toNode)
							continue;

						double edgeVal = 0.0;

						if(edgeEvals == null) {
							// If the cache processor does not have edge evaluations
							// Just check if this transition is legal based on the invalid transitions matrix
							int invalidIndex = (seqOffset + prevPos+1)*modelInfo.nPotentials;
							if(parent.invalidTransitions[invalidIndex + potential])
								continue;
						}
						else {
							// If the cache processor does have edge evaluations, then ignore the illegal transitions matrix
							// and update the expval using the edge evaluations
							FeatureEvaluation potEvals = edgeEvals[edgeIx];
							indices = potEvals.index;
							vals = potEvals.value;
							ix = 0;
							index = indices[i];
							while(index >= 0) {
								edgeVal += vals[ix] * parent.lambda[index];
								index = indices[++ix];
							}
							if(index == Short.MIN_VALUE) {
								continue;
							}
						}

						//log.debug("Stab: "+ stableBuffer.stableState[toNode]+" - " + lengthStable[toNode]);
						if(debug) {
							Assert.a(segBegin.pos == (pos-lookback-1), "Expected ",(pos-lookback-1), " was ",segBegin.pos);
						}
						double expVal = edgeVal + segBegin.mi[trans] + nodePotential;
						double prevAlpha = parent.alphas[prevPos][fromNode];
						int prevAlphaNorm = parent.alphaNorms[prevPos];

						double origBeta = segBegin.beta[fromNode];
						int origBetaNorm = segBegin.betaNorm;

						//log.debug("mi: "+buffer.mi[trans]+" s: "+nodePotential+" Base: "+(expVal - buffer.mi[trans] - nodePotential));
						int expNorm = lengthBetaHandling(seqNum, prevPos, pos, expVal, fromNode, toNode, prevAlpha, prevAlphaNorm, beta[toNode], betaNorm, nodeEvals);
						expVal -= expNorm * CleanMaximumLikelihoodSemiMarkovGradient.NORM_FACTOR;

						double transPotential = CleanMaximumLikelihoodSemiMarkovGradient.exp(expVal);
						double update = transPotential * beta[toNode];
						int updateNorm = expNorm + betaNorm;
						if(update < CleanMaximumLikelihoodSemiMarkovGradient.NORM_MIN) {
							--updateNorm;
							update *= CleanMaximumLikelihoodSemiMarkovGradient.NORM_MAX;
						}
						else if(update > CleanMaximumLikelihoodSemiMarkovGradient.NORM_MAX) {
							++updateNorm;
							update *= CleanMaximumLikelihoodSemiMarkovGradient.NORM_MIN;
						}
						// To determine the normalization value, we compare the existing beta value to the value we are about to add
						// Whichever is a larger number will dominate, and so we use that normalization value as the new value and ignore the
						// previous one.
						
						// Update the beta values
						if(updateNorm > segBegin.betaNorm) {
							//log.info(String.format("Renormalize beta from %d to %d", segBegin.betaNorm, updateNorm));
							CleanMaximumLikelihoodSemiMarkovGradient.renormalize(segBegin.beta, segBegin.betaNorm, updateNorm);
							segBegin.betaNorm = updateNorm;
						}
						else if(segBegin.betaNorm > updateNorm) {
							int expShift = updateNorm - segBegin.betaNorm;
							//log.info(String.format("Renormalize feature from %d to %d",expNorm, expNorm+expShift));
							update *= CleanMaximumLikelihoodSemiMarkovGradient.exp(expShift*CleanMaximumLikelihoodSemiMarkovGradient.NORM_FACTOR);
						}
						segBegin.beta[fromNode] += update;

						// Update expectations
						if(edgeEvals != null) {
							FeatureEvaluation potEvals = edgeEvals[edgeIx];
							indices = potEvals.index;
							vals = potEvals.value;
							ix = 0;
							index = indices[i];
							while(index != -1) {
								if(logs.expectLengthWriter != null)
									FileUtil.safeWrite(logs.expectLengthWriter, String.format("Seq %d Pos %d-%d State: %d-%d Expect #%d: %e = %e + Prob: %e * EdgeVal: %e\n", seqNum, prevPos+1, pos, modelInfo.transitionFrom[edgeIx], modelInfo.transitionTo[edgeIx], index, this.parent.expects[index]+prob*vals[i], this.parent.expects[index], prob, vals[i]));
								parent.expects[index] += prob * vals[ix];
								index = indices[++ix];
							}
						}
						
						// Updates the transition probabilities
						segBegin.transitionProb[trans] += prob;

						if(logs.betaLengthWriter != null) {
							FileUtil.safeWrite(logs.betaLengthWriter, String.format(String.format("Beta[%d][%d] = %s (%g, %d)= %s (%g, %d) + %s (%g, %d) beta[%d][%d] * %s (%g, %d) exp(Edge: %f Node: %f Stable: %f Trans: %f)\n", 
									prevPos, fromNode, CleanMaximumLikelihoodSemiMarkovGradient.printNorm(segBegin.beta[fromNode], segBegin.betaNorm), segBegin.beta[fromNode], segBegin.betaNorm, CleanMaximumLikelihoodSemiMarkovGradient.printNorm(origBeta, origBetaNorm), 
									origBeta, origBetaNorm, CleanMaximumLikelihoodSemiMarkovGradient.printNorm(beta[toNode], betaNorm), beta[toNode], betaNorm, 
									pos, toNode, CleanMaximumLikelihoodSemiMarkovGradient.printNorm(transPotential, expNorm), transPotential, expNorm, 
									edgeVal, nodePotential - stableValue, stableValue, segBegin.mi[trans]))); 
						}
					}

					// Check and renormalize the beta array.
					segBegin.betaNorm += CleanMaximumLikelihoodSemiMarkovGradient.normalize(segBegin.beta);
				}

				++lbArrayIndex;
				lengthEval = lookbacksForState[lbArrayIndex];
				lookback = lengthEval.lookback;
			}
		}
	}

	int lengthBetaHandling(int seqNum, int prevPos, int pos, double expVal, int fromNode, int toNode, double prevAlpha, int prevAlphaNorm, double betaVal, int betaNorm, FeatureEvaluation nodeEvals) {
		int norm = ((int) expVal) / CleanMaximumLikelihoodSemiMarkovGradient.NORM_FACTOR;
		expVal -= norm * CleanMaximumLikelihoodSemiMarkovGradient.NORM_FACTOR;
		
		// In addition to updating the beta array, we need to calculate a probability for this segment so we can
		// correctly calculate feature expectations
		//log.info(String.format("PrevAlpha %f Beta: %f Exp: %f zInv: %f", prevAlpha, beta[toNode], expVal, zInv));
		double afterExp = CleanMaximumLikelihoodSemiMarkovGradient.exp(expVal + CleanMaximumLikelihoodSemiMarkovGradient.NORM_FACTOR * (prevAlphaNorm + norm + betaNorm - this.parent.zNorm));
		// This if statement is here in case alpha or beta is 0, but the normalization is large.  This could cause the exp value to go to infinity and result in a NaN probability instead of 0.
		if(prevAlpha == 0 || betaVal == 0)
			prob = 0.0;
		else {
			prob = prevAlpha * betaVal * afterExp * parent.zInv;
		}
		if(Double.isNaN(prob)  || Double.isInfinite(prob)) {
			CleanMaximumLikelihoodSemiMarkovGradient.log.info(String.format("NaN = Alpha: %s * Beta: %s * Seg: %s / Z: %s",CleanMaximumLikelihoodSemiMarkovGradient.printNorm(prevAlpha, prevAlphaNorm), CleanMaximumLikelihoodSemiMarkovGradient.printNorm(betaVal, betaNorm), CleanMaximumLikelihoodSemiMarkovGradient.printNorm(CleanMaximumLikelihoodSemiMarkovGradient.exp(expVal), norm), CleanMaximumLikelihoodSemiMarkovGradient.printNorm(1/parent.zInv, parent.zNorm)));
			Assert.a(false, String.format("Seq: %d Pos: %d-%d: Bad prob (NaN) = Alpha: %e * Beta[%d] %e * %e exp(%f Norm  a:%d n:%d b:%d z:%d) * %e",
					seqNum, prevPos, pos, prevAlpha, toNode, betaVal, afterExp, expVal, prevAlphaNorm, norm, betaNorm, parent.zNorm, parent.zInv));
		}
		
		// Now update expectations for all node features for this edge
		short[] indices = nodeEvals.index;
		float[] vals = nodeEvals.value;
		int i = 0;
		short index = indices[i];
		while(index != -1) {
			if(logs.expectLengthWriter != null)
				FileUtil.safeWrite(logs.expectLengthWriter, String.format("Seq %d Pos %d-%d State: %d Expect #%d: %e = %e + Prob: %e * NodeVal: %e\n", seqNum, prevPos+1, pos, toNode, index, this.parent.expects[index]+prob*vals[i], this.parent.expects[index], prob, vals[i]));
			if(prob != 0.0)
				parent.expects[index] += prob * vals[i];
			index = indices[++i];
		}
		
		if(logs.nodeMarginalWriter != null) {
			FileUtil.safeWrite(logs.nodeMarginalWriter, String.format(
					"NodeMarg[%d][%d] = %f = %f + Alpha[%d][%d]: %s (%g n: %d) * Beta[%d][%d]: %s (%g n: %d) * seg: %s (%g n: %d) / Z: %s\n", 
					pos, toNode, nodeProb[toNode]+prob, nodeProb[toNode], prevPos, fromNode, CleanMaximumLikelihoodSemiMarkovGradient.printNorm(prevAlpha, prevAlphaNorm), prevAlpha, prevAlphaNorm, 
					pos, toNode, CleanMaximumLikelihoodSemiMarkovGradient.printNorm(betaVal, betaNorm), betaVal, betaNorm,
					CleanMaximumLikelihoodSemiMarkovGradient.printNorm(expVal, norm), expVal, norm, CleanMaximumLikelihoodSemiMarkovGradient.printNorm(1/this.parent.zInv, this.parent.zNorm)));
		}
		nodeProb[toNode] += prob;
		
		return norm;
	}
	
	void updateExpectations(int seqNum, int pos, LookbackBuffer nextLookback) {
		parent.cacheProcessor.evaluatePosition(seqNum, pos);

		// First compute the expectations for the length-dependent states.
		int invalidIndex = (seqOffset+pos)* modelInfo.nPotentials;

		for (short potential : modelInfo.orderedPotentials) {
			boolean invalid = parent.invalidTransitions[invalidIndex + potential];
			
			if (!invalid) {
				double prob = Double.NaN;
				if(potential < modelInfo.nStates) {
					prob = nodeProb[potential];
				}
				else {
					if(pos == 0)
						continue;
					prob = nextLookback.transitionProb[potential-modelInfo.nStates]; 
				}
				
				// Iterate through features for this potential.
				FeatureEvaluation potEvals = parent.evals[potential];
				short[] indices = potEvals.index;
				float[] vals = potEvals.value;
				int i = 0;
				short index = indices[i];
				while(index != -1) {
					if(logs.expectWriter != null)
						FileUtil.safeWrite(logs.expectWriter, String.format("Seq %d Pos %d Expect #%d: %e = %e + Prob: %e * Val: %e\n", seqNum, pos, index, this.parent.expects[index]+prob*vals[i], this.parent.expects[index], prob, vals[i]));
					parent.expects[index] += prob*vals[i];
					index = indices[++i];
				}
			}
		}
	}

	final void verifyMarginals(final int seqNum, final int pos, LookbackBuffer nextLookback) {	
		// As a check, we verify that the node marginals sum to one for each position.
		double sum =0.0;
		for(double x : nodeProb) {
			if(x > 1.0+CleanMaximumLikelihoodSemiMarkovGradient.ASSERTION_TOLERANCE || x < -CleanMaximumLikelihoodSemiMarkovGradient.ASSERTION_TOLERANCE)
				Assert.a(false, "Iter ",parent.iter," Seq: ",seqNum," Pos: "+pos+" Node marginals not valid "+x);
			sum += x;
		}
		if(Math.abs(1.0-sum) > CleanMaximumLikelihoodSemiMarkovGradient.ASSERTION_TOLERANCE) {
			Assert.a(false, "Iter ",parent.iter," Pos: "+pos+" Node marginals sum to "+sum+" not 1: ", ColtUtil.format(nodeProb), " at ", seqNum," ",pos);
		}

		// Verify that edge marginals sum to 1.
		if(debug && nextLookback != null) {
			double[] transitionProb = nextLookback.transitionProb;
			sum =0.0;
			for(double x : transitionProb) {
				if(x > 1+CleanMaximumLikelihoodSemiMarkovGradient.ASSERTION_TOLERANCE || x < -CleanMaximumLikelihoodSemiMarkovGradient.ASSERTION_TOLERANCE)
					Assert.a(false, "Iter ",parent.iter," Pos: "+pos+" Self-trans marginal not valid "+x);
				sum += x;
			}
	
			if(Math.abs(1.0-sum) > 0.001) {
				Assert.a(false, "Seq: ",seqNum," pos: ",pos," Edge marginals don't sum to 1.  Sum to: ", sum, ".  They are: ", ColtUtil.format(transitionProb));
			}
		}
	}	
}