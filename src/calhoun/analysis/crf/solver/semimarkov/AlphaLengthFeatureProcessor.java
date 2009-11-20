/**
 * 
 */
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
import calhoun.util.FileUtil;

final class AlphaLengthFeatureProcessor {
	static final Log log = LogFactory.getLog(CleanMaximumLikelihoodSemiMarkovGradient.class);
	static final boolean debug = log.isDebugEnabled();

	private final CleanMaximumLikelihoodSemiMarkovGradient parent;

	int seqOffset;
	int pos;
	double[] alpha;
	int alphaNorm;
	double[] stableState;
	SolverSetup modelInfo;
	final LogFiles logs;
	
	public AlphaLengthFeatureProcessor(CleanMaximumLikelihoodSemiMarkovGradient parent) {
		this.parent = parent;
		this.modelInfo = parent.modelInfo;
		this.logs = parent.logs;
	}
	
	/**
	 * In the forward pass we compute alpha values and expections. This is simpler than the backwards pass because
	 * the cache is set up for us so that we can always look at one position at a time. We have to cache previous
	 * values but we never have to look ahead.
	 */
	final void computeAlpha(final int seqNum, final int len) {
		Arrays.fill(parent.alphaNorms, Integer.MIN_VALUE);
		Arrays.fill(parent.starterAlpha, 0.0);
		double[] prevAlpha = null;
		
		seqOffset = modelInfo.seqOffsets[seqNum];

		for(pos = 0; pos < len; ++pos) {
			prevAlpha = alpha;
			alpha = parent.alphas[pos];
			alphaNorm = parent.alphaNorms[pos];
			Arrays.fill(alpha, 0.0);
			boolean alphaUpdated = false;
			if (pos == 0) {
				alphaNorm = 0;
				calcStartAlpha(alpha, seqNum);

				// Put an empty entry in the lookback so the first base has 0's initialized.
				Arrays.fill(parent.nextBuffer.stableState, 0.0);
			} else {
				parent.cacheMi(seqNum, parent.nextBuffer.mi, stableState, parent.nextBuffer.stableState, pos);
				alphaUpdated = regularAlphaUpdate(pos, parent.nextBuffer.mi, prevAlpha, alpha);
			}
			
			// Add the lookback into the array
			stableState = parent.nextBuffer.stableState;
			parent.nextBuffer = parent.lookbackBuffer.addFirst(parent.nextBuffer);
				
			// regularAlphaUpdate looks at the previous alphaNorm.
			// So, if there was an update, we need to make sure the alphaNorms agree.
			if (alphaUpdated ) {
				if (pos != 0) { 
					alphaNorm = parent.alphaNorms[pos-1]; 
				}
			}
			int norm = CleanMaximumLikelihoodSemiMarkovGradient.normalize(alpha);
			if (debug) {
				if(norm != 0) {
					if ((alphaNorm + norm > alphaNorm && norm < 0) || (alphaNorm + norm < alphaNorm && norm > 0)) {
						Assert.a(false, "Wraparound, pos=" + pos + ", norm=" + norm);
					}
				}
			}
			
			alphaNorm += norm; 

			// Now we need to loop through for the length dependent cache
			lengthAlpha(seqNum, pos);
			
			parent.alphaNorms[pos] = alphaNorm;
		}
	}
	/**
	 * Updates the alpha vector for non-length dependent states. We don't have to worry about normalization here
	 * because the regular alpha update is done before the length dependent, so these are the first values that will
	 * be set.
	 */
	private final boolean regularAlphaUpdate(final int pos, final double[] mi, final double[] lastAlpha, final double[] newAlpha) {
		double nodeVal = 0.0;
		int lastState = -1;
		boolean lengthNode = false;
		boolean ret = false;
		for (short potential : modelInfo.orderedPotentials) {
			if (potential < modelInfo.nStates) {
				if (lastState != -1) {
					if (Math.abs(nodeVal - newAlpha[lastState]) > 0.0000000000000001) {
						ret = true;
					}
					newAlpha[lastState] = nodeVal;
				}
				lastState = potential;
				nodeVal = 0.0;
				lengthNode = modelInfo.maxStateLengths[potential] > 1;
			} else {
				if (!lengthNode) {
					ret = true;
					int trans = potential - modelInfo.nStates;
					double transVal = mi[trans];
					if(!Double.isInfinite(transVal)) {
						int from = modelInfo.transitionFrom[trans];
						if(logs.alphaWriter != null)
							FileUtil.safeWrite(logs.alphaWriter, String.format("alpha[%d][%d] = %s = %s + alpha[%d][%d] %s * %s exp(%f)\n", pos, lastState, CleanMaximumLikelihoodSemiMarkovGradient.printNorm(nodeVal + lastAlpha[from] * CleanMaximumLikelihoodSemiMarkovGradient.exp(mi[trans]), alphaNorm), CleanMaximumLikelihoodSemiMarkovGradient.printNorm(nodeVal, alphaNorm), pos-1, from, CleanMaximumLikelihoodSemiMarkovGradient.printNorm(lastAlpha[from], alphaNorm), CleanMaximumLikelihoodSemiMarkovGradient.printNorm(CleanMaximumLikelihoodSemiMarkovGradient.exp(mi[trans]), 0), mi[trans]));
						nodeVal += lastAlpha[from] * CleanMaximumLikelihoodSemiMarkovGradient.exp(transVal);
					}
				}
			}
		}
		if (Math.abs(nodeVal - newAlpha[lastState]) > 0.0000000000000001) {
			ret = true;
		}
		newAlpha[lastState] = nodeVal;
		return ret;
	}

	/** Updates an alpha entry with a weighted sum of features values for a given potential */
	private final void lengthAlpha(final int seqNum, final int pos) {
		parent.cacheProcessor.evaluateSegmentsEndingAt(seqNum, pos);
		/*
		 * Updates an existing alpha by adding in: potentialValue - The value of any length-dependent features for
		 * this node f(y, i, d) and edge f(y', y, i, d) stableValue - The value of the non-length dependent node
		 * features summed across the length of this segment mis - The value of the non-length dependent transition
		 * from the previous node to this one f(y', y, i-d)
		 */
		for(int i=0; i<parent.nSemiMarkovStates; ++i) {
			LengthFeatureEvaluation[] lookbacksForState = parent.lengthEvals[i];
			StatePotentials statePotentials = modelInfo.statesWithLookback[i];
			byte toNode = statePotentials.state;
			
			int lbIndex=0;
			LengthFeatureEvaluation lengthEval = lookbacksForState[lbIndex];
			int lookback = lengthEval.lookback;
			while(lookback != -1) {
				//log.info("Pos: "+pos+"\t State: "+modelInfo.statesWithLookback[i].state+"\t Lookback: "+lookback);
				int prevPos = pos - lookback - 1;
				// For speed I hand inline RecyclingBuffer.get
				LookbackBuffer buffer = parent.lookbackBuffer.array[(parent.lookbackBuffer.currentStart+lookback)%parent.lookbackBuffer.length];

				// Handle evaluation of the node potentials
				FeatureEvaluation nodeEvals = lengthEval.nodeEval;
				short[] indices = nodeEvals.index;
				float[] vals = nodeEvals.value;
				int ix = 0;
				short index = indices[ix];
				double stableValue = stableState[toNode] - buffer.stableState[toNode];
				double nodePotential = stableValue;
				while(index >= 0) {
					nodePotential += vals[ix] * parent.lambda[index];
					index = indices[++ix];
				}
				if(debug)
					Assert.a(index != Short.MIN_VALUE, "Node lengths should only be returned in the cache if they are valid");

				if(prevPos < 0) {
					// If this is the first segment, then we don't worry about edges and handle the node directly.
					double nodeVal = nodePotential + parent.starterAlpha[toNode];
					int norm = ((int) nodeVal) / CleanMaximumLikelihoodSemiMarkovGradient.NORM_FACTOR;
					nodeVal -= norm * CleanMaximumLikelihoodSemiMarkovGradient.NORM_FACTOR;
		
					if (norm > alphaNorm) {
						CleanMaximumLikelihoodSemiMarkovGradient.renormalize(alpha, alphaNorm, norm);
						alphaNorm = norm;
						//log.info("Renormalized alpha: "+ColtUtil.format(alpha));
					} else if (norm < alphaNorm) {
						nodeVal += CleanMaximumLikelihoodSemiMarkovGradient.NORM_FACTOR * (norm - alphaNorm);
					}

					if(logs.alphaLengthWriter != null) {
						FileUtil.safeWrite(logs.alphaLengthWriter, String.format("seq: %d alpha[%d][%d] = %s = %s + %s (Pot: %f Starter: %f)\n", seqNum, pos, toNode, CleanMaximumLikelihoodSemiMarkovGradient.printNorm(alpha[toNode] + CleanMaximumLikelihoodSemiMarkovGradient.exp(nodeVal), alphaNorm), CleanMaximumLikelihoodSemiMarkovGradient.printNorm(alpha[toNode], alphaNorm), CleanMaximumLikelihoodSemiMarkovGradient.printNorm(CleanMaximumLikelihoodSemiMarkovGradient.exp(nodeVal), alphaNorm), nodePotential, parent.starterAlpha[toNode])); 
					}
					alpha[toNode] += CleanMaximumLikelihoodSemiMarkovGradient.exp(nodeVal);
				}
				else {
					// If this is not the first segment, we need to deal with edges coming into this segment
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
							int invalidIndex = (seqOffset+prevPos+1)*modelInfo.nPotentials;
							if(parent.invalidTransitions[invalidIndex + potential]) {
								//log.info("Illegal transition: "+fromNode+"-"+toNode+" at pos: "+prevPos);
								continue;
							}
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
						
						// Renormalize and update the exp value.
						double expVal = edgeVal + buffer.mi[trans] + nodePotential;
						int expNorm = ((int) expVal)/CleanMaximumLikelihoodSemiMarkovGradient.NORM_FACTOR;
						expVal -= expNorm*CleanMaximumLikelihoodSemiMarkovGradient.NORM_FACTOR;
						int prevNorm = parent.alphaNorms[prevPos];
						double prevAlpha = parent.alphas[prevPos][fromNode];
						if(prevNorm == Integer.MIN_VALUE) {
							// Previous segment begins at an invalid position and state.
							continue;
						}

						int updateNorm = expNorm + prevNorm;
						double update = CleanMaximumLikelihoodSemiMarkovGradient.exp(expVal) * prevAlpha;
						if(update == 0.0)
							continue;
						if(update < CleanMaximumLikelihoodSemiMarkovGradient.NORM_MIN) {
							--updateNorm;
							update *= CleanMaximumLikelihoodSemiMarkovGradient.NORM_MAX;
						}
						else if(update > CleanMaximumLikelihoodSemiMarkovGradient.NORM_MAX) {
							++updateNorm;
							update *= CleanMaximumLikelihoodSemiMarkovGradient.NORM_MIN;
						}
						
						if(updateNorm > alphaNorm) {
							// Our updated value is larger than the existing alpha value, renormalize that alpha vector.
							CleanMaximumLikelihoodSemiMarkovGradient.renormalize(alpha, alphaNorm, updateNorm);
							//log.info("Renormalized alpha from "+alphaNorm+" to "+updateNorm +" : "+ColtUtil.format(alpha));
							alphaNorm = updateNorm;
						}
						else if(alphaNorm > updateNorm) {
							int expShift = updateNorm - alphaNorm;
							update *= CleanMaximumLikelihoodSemiMarkovGradient.exp(expShift*CleanMaximumLikelihoodSemiMarkovGradient.NORM_FACTOR);
						}
						
						if(logs.alphaLengthWriter != null) {
							FileUtil.safeWrite(logs.alphaLengthWriter, String.format("seq: %d alpha[%d][%d] %s = %s (%g, %d) + %s (%g, %d) alpha[%d][%d] * %s (%g, %d) exp(EdgeLength: %f NodeLength: %f Edge: %f Node: %f )\n", 
									seqNum, pos, toNode, CleanMaximumLikelihoodSemiMarkovGradient.printNorm(alpha[toNode] + update, alphaNorm), CleanMaximumLikelihoodSemiMarkovGradient.printNorm(alpha[toNode], alphaNorm), alpha[toNode], alphaNorm, 
									CleanMaximumLikelihoodSemiMarkovGradient.printNorm(prevAlpha, parent.alphaNorms[prevPos]), prevAlpha, parent.alphaNorms[prevPos],
									prevPos, modelInfo.transitionFrom[trans], CleanMaximumLikelihoodSemiMarkovGradient.printNorm(CleanMaximumLikelihoodSemiMarkovGradient.exp(expVal), expNorm), CleanMaximumLikelihoodSemiMarkovGradient.exp(expVal), expNorm, edgeVal, nodePotential - stableValue, buffer.mi[trans], stableValue));
						}
									
						alpha[toNode] += update;
					}
				}
				
				++lbIndex;
				lengthEval = lookbacksForState[lbIndex];
				lookback = lengthEval.lookback;
			}
		}
	}

	/**
	 * A specialized version of calcMi for the first position in a sequence. Has the special property that constant
	 * edge features are not included.
	 */
	void calcStartAlpha(double[] currentAlpha, int seq) {
		parent.cacheProcessor.evaluatePosition(seq, 0);
		int invalidIndex = seqOffset*modelInfo.nPotentials;
		for(short potential : modelInfo.orderedPotentials) {
			if(potential < modelInfo.nStates) {
				boolean invalid = parent.invalidTransitions[invalidIndex + potential];
				double features = invalid ? Double.NEGATIVE_INFINITY : 0.0;

				// Add up all features for this potential.
				FeatureEvaluation potEvals = parent.evals[potential];
				short[] indices = potEvals.index;
				float[] vals = potEvals.value;
				int i = 0;
				short index = indices[i];
				while(index != -1) {
					if(index == Short.MIN_VALUE) {
						features = Double.NEGATIVE_INFINITY;
						break;
					}
					features += vals[i]*parent.lambda[index];
					index = indices[++i];
				}
				if(modelInfo.maxStateLengths[potential]> 1) {
					parent.starterAlpha[potential] = features;
				}
				else {
					currentAlpha[potential] = CleanMaximumLikelihoodSemiMarkovGradient.exp(features);
				}
			}
		}		
	}
}