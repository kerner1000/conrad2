/**
 * 
 */
package calhoun.analysis.crf.solver.check;

import java.io.BufferedWriter;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.CRFObjectiveFunctionGradient;
import calhoun.analysis.crf.LocalPathSimilarityScore;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.analysis.crf.scoring.SimScoreMaxStateAgreement;
import calhoun.analysis.crf.solver.CacheProcessor;
import calhoun.util.Assert;
import calhoun.util.ColtUtil;
import calhoun.util.DenseIntMatrix2D;
import calhoun.util.FileUtil;

public class CachedAOFGradient implements CRFObjectiveFunctionGradient {
	private static final Log log = LogFactory.getLog(CachedAOFGradient.class);
	private static final boolean debug = log.isDebugEnabled();

	boolean printALot = false;

	CacheProcessor cacheProcessor;
	
	LocalPathSimilarityScore score = new SimScoreMaxStateAgreement();
	
	String scoreAlphaFile = null;
	String expectedProductFile = null;
	BufferedWriter scoreAlphaWriter = null;
	BufferedWriter expectedProductWriter = null;
	
	/// Cache feature information
	short[] id;
	byte[] potentialIx;
	float[] val;

	/// Cached value of the Mi matrix for all of the features present at every position
	/// Mi is stored as a sparse matrix
	int miLength;
	double[] constMi;
	double[] mi;
	short[] transitionFrom;
	short[] transitionTo;
	short[] orderedPotentials;
	boolean[] invalidTransitions;
	int totalPositions;

	/// Cached values of the sums of each feature value through the whole training set.
	double[] featureSums;

	/// Index into the feature arrays of the first feature for each postion of each sequence.
	int[] starts;

	/// Index into the starts array of the first position of each sequence.
	int[] seqOffsets;

	/// Number of sequences in the training data set
	int nSeqs; 
	
	/// Number of constant features.
	int nConstantFeatures; 
	
	List<? extends TrainingSequence<?>> data;
	ModelManager fm;
	int nFeatures;
	int nStates;
	int nPotentials;
	int nTransitions;
	int iter = 0;

	double[] prevAlpha;
	double[] alpha;
	double[][] betas;
	double[] betaNorms;
	double[] expects;
	//double[][] allAlpha;
	double[][] scoreAlpha;
	double[][] scoreBeta;
	//double[] alphaNorms;
	double[][] edgeProb;
	double[][] nodeProb;
	private DenseIntMatrix2D transitionIndex;
	double[] localExpects;
	double[] temp1,temp2;


	double exp(double val1) {
		return Math.exp(val1);
	}
	
	double log(double val1) {
		return Math.log(val1);
	}
	
	boolean allPaths;
	public void setAllPaths(boolean allPaths) {
		this.allPaths = allPaths;
	}
	
	public void setTrainingData(ModelManager fm, List<? extends TrainingSequence<?>> data) {
		this.fm = fm;
		this.data = data;
		nFeatures = fm.getNumFeatures();
		nStates = fm.getNumStates();
		nSeqs = data.size();

		expects = new double[nFeatures];
		localExpects = new double[nFeatures];
		prevAlpha = new double[nStates];
		alpha = new double[nStates];

		FeatureCache cache = new FeatureCache(fm, data, allPaths);
		Assert.a(nStates == cache.nStates);
		
		// Create local references to cache values
		orderedPotentials = cache.orderedPotentials;
		id = cache.id;
		potentialIx = cache.potentialIx;
		val = cache.val;

		transitionFrom = cache.transitionFrom;
		transitionTo = cache.transitionTo;
		transitionIndex = cache.transitionIndex;
		miLength = cache.nTransitions;
		constMi = new double[miLength];
		mi = new double[miLength];

		featureSums = cache.featureSums;
		starts = cache.starts;
		seqOffsets = cache.seqOffsets;
		invalidTransitions = cache.invalidTransitions;
		nConstantFeatures = cache.numConstantFeatures;
		nPotentials = cache.nPotentials;
		nTransitions = cache.nTransitions;
		totalPositions = cache.totalPositions;
		
		// Initialize betas (for use later, in the gradient computation)
		betas = new double[cache.longestSeq][nStates];
		betaNorms = new double[cache.longestSeq];
		
		// Initialize variables for use in computing gradient of AOF:
		//allAlpha = new double[cache.longestSeq][nStates];
		//alphaNorms = new double[cache.longestSeq];
		scoreAlpha = new double[cache.longestSeq][nStates];
		scoreBeta = new double[cache.longestSeq][nStates];
		edgeProb = new double[cache.longestSeq][cache.nTransitions]; //[pos][yprev][y], not defined for pos=0
		nodeProb = new double[cache.longestSeq][nStates];
		temp1 = new double[nStates];
		temp2 = new double[nStates];		
	}
	
	public void clean() {
	}
	
	/* Returns the value of the objective function E_w(S(y,y0,x)), and is passed a vector grad which is populated with
	 * the gradient of this objective function.
	 */
	public double apply(double[] param, double[] grad) {
		scoreAlphaWriter = FileUtil.safeOpen(scoreAlphaFile);
		expectedProductWriter = FileUtil.safeOpen(expectedProductFile);

		Arrays.fill(grad, 0);
		double result = 0.0;
		double[] seqGrad = new double[grad.length];
		
		for (int jfeat=0; jfeat<nFeatures; jfeat++) {
			Assert.a(!Double.isNaN(grad[jfeat]));
		}		
		
		// Calculate the constant Mi matrix
		// CalcMi produces an exponentiated matrix, so we take the log the get the right init values
		// and copy them to the constMi matrix.
		Arrays.fill(constMi, 0.0);
		calcMi(-1, 0, starts[0], param);
		for(int i = 0; i<miLength; ++i) {
			constMi[i] = log(mi[i]);
		}
		
		// Iterate through sequences
		Arrays.fill(expects, 0);
		int seqStart = 0;
		for(int i = 0; i<nSeqs; ++i) {
			Arrays.fill(seqGrad, 0);
			int len = seqOffsets[i+1] - seqOffsets[i]; 
			
			// Work backwards, computing betas.
			Arrays.fill(betas[len-1], 1.0);
			betaNorms[len-1] = 0;
			int cacheStop = starts[seqStart + len];
			for (int pos = len - 1; pos > 0; --pos) {
				int overallPosition = seqStart + pos;
				int cacheStart = starts[overallPosition];
				calcMi(overallPosition, cacheStart, cacheStop, param);
				cacheStop = cacheStart;
				quickBetaUpdate(betas[pos], betas[pos-1]);
				
				//mi.zMult(betas[j], betas[j-1], 1, 0, false);
				double n = normalizePotential(betas[pos-1]);
				betaNorms[pos-1] = betaNorms[pos] + log(n);
				Assert.a(!Double.isNaN(log(n)));
			}
			// after above, betaNorms and beta are ready for use; exp(betaNorms)*betas is what appears in my
			// mathematical formulae as "beta".
			
			
			// Now work forwards
			double logZ = Double.NEGATIVE_INFINITY;	// This should always get initialized.  Blow up if not.
			double alphaNorm = 0;
			double prevAlphaNorm = 0;
			int cacheStart = starts[seqStart];
			for (int j=0; j<nFeatures; j++) { localExpects[j] = 0.0; }
			for(int pos=0; pos<len; ++pos) {
				int overallPosition = seqStart + pos;
				double[] beta = betas[pos];
				double betaNorm = betaNorms[pos];
				cacheStop = starts[overallPosition+1];
				if(pos == 0) {
					calcStartAlpha(overallPosition, cacheStart, cacheStop, param);
					alphaNorm = log(normalizePotential(alpha));
					
					// We now have everything needed to compute Z.
					logZ = log(ColtUtil.dotProduct(alpha, beta)) + betaNorm + alphaNorm;
				}
				else {
					calcMi(overallPosition, cacheStart, cacheStop, param);
					//r * M (or M'*r)
					//mi.zMult(prevAlpha, alpha, 1, 0, true);
					quickAlphaUpdate(prevAlpha, alpha);
					alphaNorm = prevAlphaNorm + log(normalizePotential(alpha));
					
					// Verify our calculations by checking the Z.
					// Expensive, so only enable as needed during debugging
					//double newZ = log(ColtUtil.dotProduct(alpha, beta)) + betaNorm + alphaNorm;
					//Assert.a(Math.abs(newZ-logZ) < 0.0000001*Math.abs(logZ), "New Z:",newZ," Old was: ", logZ);
				}
				
				
				// Iterate through the potentials to update feature expectations
				double nodeNorm = exp(alphaNorm + betaNorm - logZ);
				if (Double.isNaN(nodeNorm)) {
					Assert.a(!Double.isNaN(nodeNorm)," alphaNorm = " + alphaNorm + " betaNorm " + betaNorm + "  logZ " + logZ);
				}
				double edgeNorm = exp(prevAlphaNorm + betaNorm - logZ);
				if (Double.isNaN(nodeNorm) || Double.isNaN(edgeNorm) || Double.isInfinite(nodeNorm) || Double.isInfinite(edgeNorm) ) {
					Assert.a(false,"nodeNorm = " + nodeNorm + "  edeNorm = " + edgeNorm + "  alphaNorm = " + alphaNorm + "  betaNorm = " + betaNorm + "  prevAlphaNorm = " + prevAlphaNorm + "  logZ = " + logZ);
				}
				updateExpectations(overallPosition, pos, pos != 0, cacheStart, cacheStop, nodeNorm, edgeNorm, beta);

				if(debug) {
					if((i < 2 || i == nSeqs-1) && (pos < 2 || pos >= len-2)) {
						log.debug(String.format("Pos: %d expects: %s alphas: %s (norm %f) betas: %s (norm %f)", pos, ColtUtil.format(expects), ColtUtil.format(alpha), alphaNorm, ColtUtil.format(beta), betaNorm));
					}
				}

				// Recycle the arrays
				double[] swap = prevAlpha;
				prevAlpha = alpha;
				alpha = swap;
				prevAlphaNorm = alphaNorm;
				cacheStart = cacheStop;
			}
			// after above, alphaNorms and allAlpha are ready for use; exp(alphaNorms)*allAlpha is what appears in my
			// mathematical formulae as "alpha".

			// Mathematically, it should be the case that
			// 1) sum_yprev edgeprob[pos][yprev][y] = nodeprob[pos][y]
			// 2) sum_y edgeprob[pos][yprev][y] = nodeprob[pos-1][yprev]
			// But this may fail due to round-off errors.  In fact, the nodeprob could be zero
			// even as some of the edgeprobs are nonzero, perhaps leading to a division by zero
			// and NaN troubles that snowball.
			for (int pos = 1; pos <len; pos++) {
				Arrays.fill(temp1,0);
				Arrays.fill(temp2,0);
				for (int trans=0; trans<nTransitions; trans++) {
					int yprev = transitionFrom[trans];
					int y = transitionTo[trans];
					//double eps2 = 0.00000000001;
					double eps2 = Double.MIN_VALUE;
					if ((edgeProb[pos][trans]<eps2) || (Double.isNaN(edgeProb[pos][trans]))) { edgeProb[pos][trans] = eps2; }
					double ep = edgeProb[pos][trans];
					if (!(ep>=0)) {
						Assert.a(false,"ep = " + ep + "   eps = " + eps2);
					}
					temp1[yprev] += ep;
					temp2[y] += ep;
				}
				for (int j=0; j<nStates; j++) {			
					if (temp1[j] > nodeProb[pos-1][j]) { nodeProb[pos-1][j] = temp1[j];  }
					if (temp2[j] > nodeProb[pos][j]) { nodeProb[pos][j] = temp2[j]; }
				}
				for (int trans=0; trans<nTransitions; trans++) {
					int yprev = transitionFrom[trans];
					int y = transitionTo[trans];
					double ep = edgeProb[pos][trans];
					if (nodeProb[pos-1][yprev] < ep) {
						Assert.a( false , " ep = " + ep + " np = " + nodeProb[pos-1][yprev]);
					}
					Assert.a( nodeProb[pos][y] >= ep );
				}
			}
			
			for (int pos = 0; pos <len; pos++) {
				for (int stat = 0; stat<nStates; stat++) {
					double np = nodeProb[pos][stat];
					if (printALot) { System.out.println("At pos=" + pos + "  the node to y=" + stat + " has probability " + np);} 					
				}
				for (int trans=0; trans<nTransitions; trans++) {
					int yprev = transitionFrom[trans];
					int y = transitionTo[trans];
					double ep = edgeProb[pos][trans];
					if (printALot) { System.out.println("At pos=" + pos + "  the edge from yprev="+yprev + " to y=" + y + " has probability " + ep); }
				}
			}
			
			
			// Now do another backward pass
			// scoreBeta is being defined from 0 to len-1
			Arrays.fill(scoreBeta[len-1], 0.0);
			for (int pos = len - 1; pos > 0; --pos) {
				for (int yprev=0; yprev<nStates; yprev++) {
					scoreBeta[pos-1][yprev]=0.0;
				}
				for (int trans=0; trans<nTransitions; trans++) {
					int yprev = transitionFrom[trans];
					int y = transitionTo[trans];
					scoreBeta[pos-1][yprev] += edgeProb[pos][trans]*score.evaluate(yprev,y,data.get(i),pos);
					double ep = edgeProb[pos][trans];
					if (ep>0) {
						double np = nodeProb[pos][y]; 
						Assert.a(np>0);
						double ratio = ep / np;
						if (!Double.isNaN(ratio)) {
							scoreBeta[pos-1][yprev] += (edgeProb[pos][trans]/nodeProb[pos][y])*scoreBeta[pos][y];
						} else {
							Assert.a(false, "   ep = " + ep + "   np = " + np);
						}
					}
				}
			}
			
			if (printALot) {
				// Now print out lots of diagnostic info
				for (int pos=0; pos<len; pos++ ) {
					for (int y=0; y<nStates; y++) {
						System.out.println("scoreAlpha at position " + pos + " and state " + y + " is " + scoreAlpha[pos][y]);
						System.out.println("scoreBeta at position " + pos + " and state " + y + " is " + scoreBeta[pos][y]);
					}
				}
			}
			
			
			// Now do another forward pass
			// scoreAlpha is being defined from 0 to len-1
			Arrays.fill(scoreAlpha[0],0.0); // should actually be alpha_0(y_0) = nodeprob(0,y0)*score(start,y0,start,y00,x,0), but I don't know how to do the start just yet.
			for (int pos=1; pos<len; pos++) {
				for (int y=0; y<nStates; y++) {
					scoreAlpha[pos][y] = 0.0;
				}
				for (int trans=0; trans<nTransitions; trans++) {
					double ep = edgeProb[pos][trans];
					if (ep>0) {
						int yprev = transitionFrom[trans];
						int y = transitionTo[trans];
						double update = (ep/nodeProb[pos-1][yprev])*scoreAlpha[pos-1][yprev];
						update += ep*score.evaluate(yprev,y,data.get(i),pos);
						scoreAlpha[pos][y] += update;
	
						if(scoreAlphaWriter != null) {
							FileUtil.safeWrite(scoreAlphaWriter, String.format("Seq: %d alpha[%d][%d] = %g = %g + Pr: %g * alpha[%d][%d] %g + Pr: %g * Score: %g\n",
									i, pos, y, scoreAlpha[pos][y], scoreAlpha[pos][y]-update, ep/nodeProb[pos-1][yprev], pos-1, yprev, scoreAlpha[pos-1][yprev], ep, score.evaluate(yprev,y,data.get(i),pos)));
						}
					}
				}	
			}
			
			for (int jfeat=0; jfeat<nFeatures; jfeat++) {
				Assert.a(!Double.isNaN(grad[jfeat]));
			}		
			
			// Now finish up; calculate the objective function and its gradient.
			double thisresult = 0.0;
			for (int pos=0; pos<len; pos++) {
				if (pos>0) {
					for (int trans=0; trans<nTransitions; trans++) {
						int yprev = transitionFrom[trans];
						int y = transitionTo[trans];
						thisresult += edgeProb[pos][trans] * score.evaluate(yprev,y,data.get(i),pos);	
					}
				}
				int overallPosition = seqStart + pos;
				cacheStart = starts[overallPosition];
				cacheStop = starts[overallPosition+1];
				Assert.a(nFeatures>0);
				for (int jfeat=0; jfeat<nFeatures; jfeat++) {
					Assert.a(!Double.isNaN(grad[jfeat]));
				}	
				for (int jfeat=0; jfeat<nFeatures; jfeat++) {
					Assert.a(!Double.isNaN(grad[jfeat]));
				}	
				updateGrad(overallPosition,pos,i,cacheStart,cacheStop,seqGrad);
				if (printALot) {
					System.out.println("After adding in pos=" + pos + " the gradient is " + grad[0] + " , " + grad[1]);
				}
			}
			for (int jfeat=0; jfeat<nFeatures; jfeat++) {
				Assert.a(!Double.isNaN(grad[jfeat]));
			}
			if (printALot) {
				System.out.println("For sequence number " + i + " the expectation of S is " + thisresult);
			}
			for (int jfeat=0; jfeat<nFeatures; jfeat++) {
				if (printALot) { 
					System.out.println("For seq " + i + " and feature " + jfeat + " the local expect is " + localExpects[jfeat] + "  and expects is " + expects[jfeat]);
				}
				grad[jfeat] += seqGrad[jfeat] - localExpects[jfeat] * thisresult;
			}
			for (int jfeat=0; jfeat<nFeatures; jfeat++) {
				Assert.a(!Double.isNaN(grad[jfeat]));
			}		
			
			if(debug) {
				log.debug(String.format("Iter: %d Seq: %d Expected Score: %g Grad: %s Expected Features: %s Expected Product: %s", iter, i, thisresult, ColtUtil.format(grad), ColtUtil.format(localExpects), ColtUtil.format(seqGrad)));
			}
			result += thisresult;
			seqStart += len;
		}
		
		if(log.isInfoEnabled()) {
			// Report average per-label numbers.
			log.info(String.format("Iter: %d Val: %g Grad: %s", iter, result, ColtUtil.format(grad)));
		}

		Assert.a(!Double.isNaN(result));
		for (int jfeat=0; jfeat<nFeatures; jfeat++) {
			Assert.a(!Double.isNaN(grad[jfeat]));
		}		
		result = result/totalPositions; 
		for(int i=0; i<grad.length; ++i) {
			grad[i] = grad[i]/totalPositions;
		}
		iter += 1;
		FileUtil.safeClose(scoreAlphaWriter);
		FileUtil.safeClose(expectedProductWriter);
		return result;
	}

	void updateGrad(int overallPos, int pos, int seqI, int posCurrent, int posStop, double[] grad) {

		for (int jfeat=0; jfeat<nFeatures; jfeat++) {
			Assert.a(!Double.isNaN(grad[jfeat]));
		}	
		
		int constCurrent = 0;
		
		
		for (int j=0; j<nFeatures; j++) {
			if(Double.isNaN(grad[j])) {
				Assert.a(false, "j= " + j);
			}
		}
		
		// Constant features
		short constId = -1;
		byte constPotential = -1;
		double constVal = Double.NaN;
		if(constCurrent < nConstantFeatures) {
			constId = id[constCurrent];
			constPotential = potentialIx[constCurrent];
			constVal = val[constCurrent];
			++constCurrent;
		}
		
		// Positional features
		short posId = -1;
		byte posPotential = -1;
		double posVal = Double.NaN;
		if(posCurrent < posStop) {
			posId = id[posCurrent];
			posPotential = potentialIx[posCurrent];
			posVal = val[posCurrent];
			++posCurrent;
		}

		int invalidIndex = overallPos*nPotentials;
		for(short potential : orderedPotentials) {
			boolean invalid = invalidTransitions[invalidIndex + potential];
			
			// Include constant features for this potential.  Need to skip over constant features for invalid potentials
			while(constPotential == potential) {
				if(!invalid && ((pos>0) || potential < nStates)) {
					// DO SOMETHING RIGHT HERE
					//	expects[constId] += prob*constVal;
					// constId = index of the relevant feature
					// constVal = value of the relevant fetaure
					// constPotential = index of relevant potential; at this place it is equal to potential
					if (potential < nStates) {
						// It's a node
						int y = potential;
						if (pos>0) {
							for (int yprev=0; yprev<nStates; yprev++) {
								int trans = transitionIndex.getQuick(yprev,y);
								if (trans >= 0) {
									Assert.a(yprev == transitionFrom[trans]);
									Assert.a(y == transitionTo[trans]);							
									double inner = 0;
									double ep = edgeProb[pos][trans];
									if (ep>0) {
										inner += (edgeProb[pos][trans]/nodeProb[pos-1][yprev])*scoreAlpha[pos-1][yprev];
										inner += (edgeProb[pos][trans])*score.evaluate(yprev,y,data.get(seqI),pos);
										inner += (edgeProb[pos][trans]/nodeProb[pos][y])*scoreBeta[pos][y];
										grad[constId] += constVal*inner;
										if(expectedProductWriter != null) {
											FileUtil.safeWrite(expectedProductWriter, String.format("Seq: %d Pos: %d State: %d\tFeat: %d = %g = %g + Val: %g * (s: %g * ep: %g + a: %g * ms: %g + b: %g * me: %g)\n",
													seqI, pos, y, constId, grad[constId], grad[constId]-inner*constVal, constVal, 
													score.evaluate(yprev,y,data.get(seqI),pos), (edgeProb[pos][trans]), scoreAlpha[pos-1][yprev], (edgeProb[pos][trans]/nodeProb[pos-1][yprev]),
													scoreBeta[pos][y], edgeProb[pos][trans]/nodeProb[pos][y]));
										}
									}
									if(Double.isNaN(grad[constId])) {
										Assert.a(false, "posVal = " + posVal + "  inner = " + inner + "pos = " + pos + "   trans = " + trans);
									}
								}
							}
						} else {
							double inner = scoreBeta[pos][y];
							grad[constId] += constVal*inner;					
							if(expectedProductWriter != null) {
								FileUtil.safeWrite(expectedProductWriter, String.format("Seq: %d Pos: %d State: %d\tFeat: %d = %g = %g + Val: %g * Beta[%d][%d]: %g:\n",
										seqI, pos, y, constId, grad[constId], grad[constId]-inner*constVal, constVal, pos, y, inner));
							}
						}
					} else {
						// it's an edge
						int trans = potential - nStates;
						int yprev = transitionFrom[trans];
						int y = transitionTo[trans];
						double inner = 0;
						double ep = edgeProb[pos][trans];
						if (ep>0) {
							inner += (ep/nodeProb[pos-1][yprev])*scoreAlpha[pos-1][yprev];
							inner += (ep)*score.evaluate(yprev,y,data.get(seqI),pos);
							inner += (ep/nodeProb[pos][y])*scoreBeta[pos][y];
							grad[constId] += constVal*inner;
							if(expectedProductWriter != null) {
								FileUtil.safeWrite(expectedProductWriter, String.format("Seq: %d Pos: %d Edge: %d-%d\tFeat: %d = %g = %g + Val: %g * (s: %g * ep: %g + a: %g * ms: %g + b: %g * me: %g)\n",
										seqI, pos, yprev, y, constId, grad[constId], grad[constId]-inner*constVal, constVal,
										score.evaluate(yprev,y,data.get(seqI),pos), ep, scoreAlpha[pos-1][yprev], ep/nodeProb[pos-1][yprev],
										scoreBeta[pos][y], ep/nodeProb[pos][y]));
							}
						}
						if(Double.isNaN(grad[constId])) {
							Assert.a(false, "posVal = " + posVal + "  inner = " + inner + "pos = " + pos + "   trans = " + trans);
						}
					}
				}
				if(constCurrent < nConstantFeatures) {
					constId = id[constCurrent];
					constVal = val[constCurrent];
					constPotential = potentialIx[constCurrent];
					++constCurrent;
				}
				else { break; }
			}
			
			for (int j=0; j<nFeatures; j++) {
				if(Double.isNaN(grad[j])) {
					Assert.a(false, "j= " + j);
				}
			}
			
			// Include cached features for this potential.  The cache shoudl never have features for invalid potentials.
			if(!invalid) {
				while(posPotential == potential) {
					// DO SOMETHING RIGHT HERE
					if (potential < nStates) {
						// It's a node
						int y = potential;
						if (pos>0) {
							for (int yprev=0; yprev<nStates; yprev++) {
								int trans = transitionIndex.getQuick(yprev,y);
								if (trans >= 0) {
									Assert.a(yprev == transitionFrom[trans]);
									Assert.a(y == transitionTo[trans]);							
									double inner = 0;
									double ep = edgeProb[pos][trans];
									if (ep>0) {
										Assert.a(nodeProb[pos-1][yprev]>=(ep-0.00000001));
										Assert.a(nodeProb[pos-1][yprev]>0);
										Assert.a(nodeProb[pos][y]>=(ep-0.00000001));
										Assert.a(nodeProb[pos][y]>0);									
										inner += (ep/nodeProb[pos-1][yprev])*scoreAlpha[pos-1][yprev];
										inner += (ep)*score.evaluate(yprev,y,data.get(seqI),pos);
										inner += (ep/nodeProb[pos][y])*scoreBeta[pos][y];
										if(Double.isNaN(inner)) {
											Assert.a(false, "npyprev = " + nodeProb[pos-1][yprev] + "  npy = " + nodeProb[pos][y] + "pos = " + pos + "   ep = " + ep + "  score = " + score.evaluate(yprev,y,data.get(seqI),pos) + "  alpha= "+ scoreAlpha[pos-1][yprev] + "   beta=" + scoreBeta[pos][y]);
										}
										grad[posId] += posVal*inner;
										if(expectedProductWriter != null) {
											FileUtil.safeWrite(expectedProductWriter, String.format("Seq: %d Pos: %d State: %d\tFeat: %d = %g = %g + Val: %g * (s: %g * ep: %g + a: %g * ms: %g + b: %g * me: %g)\n",
													seqI, pos, y, posId, grad[posId], grad[posId]-inner*posVal, posVal,
													score.evaluate(yprev,y,data.get(seqI),pos), (edgeProb[pos][trans]), scoreAlpha[pos-1][yprev], (edgeProb[pos][trans]/nodeProb[pos-1][yprev]),
													scoreBeta[pos][y], edgeProb[pos][trans]/nodeProb[pos][y]));
										}
									}
									if(Double.isNaN(grad[posId])) {
										Assert.a(false, "posVal = " + posVal + "  inner = " + inner + "pos = " + pos + "   trans = " + trans);
									}
								}
							}
						} else {
							double inner = scoreBeta[pos][y];
							grad[posId] += posVal*inner;
						}
					} else {
						// it's an edge
						int trans = potential - nStates;
						int yprev = transitionFrom[trans];
						int y = transitionTo[trans];
						double inner = 0;
						double ep = edgeProb[pos][trans];
						if (ep>0) {
							inner += (ep/nodeProb[pos-1][yprev])*scoreAlpha[pos-1][yprev];
							inner += (ep)*score.evaluate(yprev,y,data.get(seqI),pos);
							inner += (ep/nodeProb[pos][y])*scoreBeta[pos][y];
							grad[posId] += posVal*inner;
							if(expectedProductWriter != null) {
								FileUtil.safeWrite(expectedProductWriter, String.format("Seq: %d Pos: %d Edge: %d-%d\tFeat: %d = %g = %g + Val: %g * (s: %g * ep: %g + a: %g * ms: %g + b: %g * me: %g)\n",
										seqI, pos, yprev, y, posId, grad[posId], grad[posId]-inner*posVal, posVal,
										score.evaluate(yprev,y,data.get(seqI),pos), (edgeProb[pos][trans]), scoreAlpha[pos-1][yprev], (edgeProb[pos][trans]/nodeProb[pos-1][yprev]),
										scoreBeta[pos][y], edgeProb[pos][trans]/nodeProb[pos][y]));
							}
						}
						if(Double.isNaN(grad[posId])) {
							Assert.a(false, "posVal = " + posVal + "  inner = " + inner + "pos = " + pos + "   trans = " + trans);
						}
					}
				//	expects[posId] += prob*posVal;
					if(posCurrent < posStop) {
						posId = id[posCurrent];
						posVal = val[posCurrent];
						posPotential = potentialIx[posCurrent];
						++posCurrent;
					}
					else { break; }
				}
			}
		}
		
		for (int jfeat=0; jfeat<nFeatures; jfeat++) {
			Assert.a(!Double.isNaN(grad[jfeat]));
		}		
		
		// verify that we actually checked all the cache entries.
		Assert.a(constCurrent == nConstantFeatures);
		Assert.a(posCurrent == posStop);
	}
	
	
	void updateExpectations(int overallPos, int localPos, boolean includeEdges, int posCurrent, int posStop, double nodeNorm, double edgeNorm, double[] beta) {
		int constCurrent = 0;
		
		// Constant features
		short constId = -1;
		byte constPotential = -1;
		double constVal = Double.NaN;
		if(constCurrent < nConstantFeatures) {
			constId = id[constCurrent];
			constPotential = potentialIx[constCurrent];
			constVal = val[constCurrent];
			++constCurrent;
		}
		
		// Positional features
		short posId = -1;
		byte posPotential = -1;
		double posVal = Double.NaN;
		if(posCurrent < posStop) {
			posId = id[posCurrent];
			posPotential = potentialIx[posCurrent];
			posVal = val[posCurrent];
			++posCurrent;
		}

		int currentNode = -1;
		double currentBeta = 0.0f;
		int invalidIndex = overallPos*nPotentials;
		for(short potential : orderedPotentials) {
			boolean invalid = invalidTransitions[invalidIndex + potential];

			double prob = 0.0;
			if(potential < nStates) {
				currentNode = potential;
				currentBeta = beta[currentNode];
				if(!invalid) {
					prob = alpha[currentNode] * currentBeta * nodeNorm;
				}
				if (Double.isNaN(prob)) {
					Assert.a(false," alpha" + alpha[currentNode] + " beta=" + currentBeta + "  nodeNorm = " + nodeNorm);
				}
				nodeProb[localPos][currentNode] = prob;
			}
			else {
				int trans = potential - nStates;
				int yprev = transitionFrom[trans];
				if(!invalid) {
					prob = prevAlpha[yprev] * mi[trans] * currentBeta * edgeNorm;
				}
				if (Double.isNaN(prob)) {
					Assert.a(false,"prob="+ prob+"  prevAlpha="+prevAlpha[yprev]+"  mi="+mi[trans]+"  currentBeta="+currentBeta+ "  edgeNorm="+edgeNorm);
				}
				edgeProb[localPos][trans] = prob;
			}
			
			// Include constant features for this potential.  Need to skip over constant features for invalid potentials
			while(constPotential == potential) {
				if(!invalid && (includeEdges || potential < nStates)) {
					expects[constId] += prob*constVal;
					localExpects[constId] += prob*constVal;
				}
				if(constCurrent < nConstantFeatures) {
					constId = id[constCurrent];
					constVal = val[constCurrent];
					constPotential = potentialIx[constCurrent];
					++constCurrent;
				}
				else { break; }
			}

			// Include cached features for this potential.  The cache shoudl never have features for invalid potentials.
			if(!invalid) {
				while(posPotential == potential) {
					expects[posId] += prob*posVal;
					localExpects[posId] += prob*posVal;					
					if(posCurrent < posStop) {
						posId = id[posCurrent];
						posVal = val[posCurrent];
						posPotential = potentialIx[posCurrent];
						++posCurrent;
					}
					else { break; }
				}
			}
		}
		// verify that we actually checked all the cache entries.
		Assert.a(constCurrent == nConstantFeatures);
		Assert.a(posCurrent == posStop);
	}
	
	/** This is one of the most time critical parts of the entire solver.  The goal is to update the transition matrix.
	 * This function makes a lot of assumptions in order to maximize performance.
	 * 
	 * To maximize performance, we want to make one pass through the Mi matrix, setting each entry to its correct value.
	 * The value for each entry is the exponent of the sum of the weighted feature values of the edge for that entry and 
	 * its corresponding node.  The entry s0,s1 consists of the s0,s1 edge and the s1 node. 
	 * 
	 * Because node features are applied to more than 1 entry in the matrix, we use a sorting of all of the features where
	 * each node preceeds all its corresponding edges.  This allows us to keep track of only 1 node value at a time and easily
	 * apply it to all its edge features.
	 * 
	 * As we evaluate each potential, we check the cache to see if it is valid at this position and to get any features values.
	 * Note that this function very much depends on the fact that the entries in the cache will be in the correct order.
	 * 
	 * The other wrinkle is that for features that always occur (constant features), we pull them from the constant mi array, not
	 * from the cache.
	 * 
	 *  This function is also used to calculate the Mi Matrix for the constant features.
	 *  */
	void calcMi(int overallPosition, int current, int stop, double[] lambda) {
		// Features will always be in the order of the potentials.  We loop through the potentials, grabbing the features for each.
		byte cachedPotential = -1;
		double cachedVal = Double.NaN;
		if(current < stop) {
			cachedPotential = potentialIx[current];
			cachedVal = val[current]*lambda[id[current]];
			++current;
		}
		double nodeVal = Double.NaN;
		int invalidIndex = overallPosition*nPotentials;
		for(short potential : orderedPotentials) {
			boolean invalid = overallPosition != -1 && invalidTransitions[invalidIndex + potential];
			double features = invalid ? Double.NEGATIVE_INFINITY : 0.0;

			// Add up all features for this potential.
			while(cachedPotential == potential) {
				features += cachedVal;
				//Assert.a(!Double.isNaN(features));
				if(current < stop) {
					cachedVal = val[current]*lambda[id[current]];
					cachedPotential = potentialIx[current];
					++current;
				}
				else { break; }
			}

			if(potential < nStates) {
				nodeVal = features;
			}
			else {
				//log.debug(String.format("Mi[%d, %d] = %f, adding in %f to get %f", feat.yprev(), feat.y(), val, feat.value(), val*exp(feat.value()*param[feat.index()])));
				int transition = potential - nStates;
				mi[transition] = exp(features + nodeVal + constMi[transition]);
			}
		}
		// verify that we actually checked all the cache entries.
		if ( current != stop ) {
			Assert.a(false, "Pos: ", overallPosition, " Expected ", stop, " features only found ", current);
		}
	}
	
	/** A specialized version of calcMi for the first position in a sequence. Has the special property that 
	 * constant edge features are not included.  Also optionally allows you to set a node value used to initialize the alphas. */
	void calcStartAlpha(int overallPosition, int posCurrent, int posStop, double[] lambda) {
		int constCurrent = 0;
		
		// Constant features
		byte constPotential = -1;
		double constVal = Double.NaN;
		if(constCurrent < nConstantFeatures) {
			constPotential = potentialIx[constCurrent];
			constVal = val[constCurrent]*lambda[id[constCurrent]];
			++constCurrent;
		}
		
		// Positional features
		byte posPotential = -1;
		double posVal = Double.NaN;
		if(posCurrent < posStop) {
			posPotential = potentialIx[posCurrent];
			posVal = val[posCurrent]*lambda[id[posCurrent]];
			++posCurrent;
		}

		int invalidIndex = overallPosition*nPotentials;
		for(short potential : orderedPotentials) {
			boolean invalid = invalidTransitions[invalidIndex + potential];
			double features = invalid ? Double.NEGATIVE_INFINITY : 0.0;

			// Find constant features for this potential.  Ignore edge features.
			while(constPotential == potential) {
				if(potential < nStates) {
					features += constVal;
					Assert.a(!Double.isNaN(features));
				}
				if(constCurrent < nConstantFeatures) {
					constVal = val[constCurrent]*lambda[id[constCurrent]];
					constPotential = potentialIx[constCurrent];
					++constCurrent;
				}
				else { break; }
			}

			if(potential < nStates) {
				// Add in cached features
				while(posPotential == potential) {
					features += posVal;
					if(posCurrent < posStop) {
						posVal = val[posCurrent]*lambda[id[posCurrent]];
						posPotential = potentialIx[posCurrent];
						++posCurrent;
					}
					else { break; }
				}
				alpha[potential] = exp(features);
			}
		}
		// verify that we actually checked all the cache entries.
		Assert.a(constCurrent == nConstantFeatures);
		Assert.a(posCurrent == posStop);
	}
	
	private double normalizePotential(double[] vec) {
		double norm = 0.0;
		int len = vec.length;
		for(int i=0; i<len; ++i) {
			norm += vec[i];
		}
		Assert.a(norm>0);
		double mult = 1/norm;
		for(int i=0; i<len; ++i) {
			vec[i] *= mult;
		}
		return norm;
	}

	private void quickBetaUpdate(double[] lastBeta, double[] newBeta) {
		Arrays.fill(newBeta, 0);
		double nodeVal = 0.0;
		for(short potential : orderedPotentials) {
			if(potential < nStates) {
				nodeVal = lastBeta[potential];
			}
			else {
				int trans = potential - nStates;
				int from = transitionFrom[trans];
				newBeta[from] += mi[trans]*nodeVal;
			}
		}
	}

	private void quickAlphaUpdate(double[] lastAlpha, double[] newAlpha) {
		double nodeVal = 0.0;
		int lastState = -1;
		for(short potential : orderedPotentials) {
			if(potential < nStates) {
				if(lastState != -1) {
					newAlpha[lastState] = nodeVal;
				}
				lastState = potential;
				nodeVal = 0.0;
			}
			else {
				int trans = potential - nStates;
				int from = transitionFrom[trans];
				nodeVal += lastAlpha[from]*mi[trans];
			}
		}
		newAlpha[lastState] = nodeVal;
	}

	public CacheProcessor getCacheProcessor() {
		return cacheProcessor;
	}

	public void setCacheProcessor(CacheProcessor cacheProcessor) {
		this.cacheProcessor = cacheProcessor;
	}

	public void setLocalPathSimilarityScore(LocalPathSimilarityScore s) {
		score = s;
	}
	
	public LocalPathSimilarityScore getLocalPathSimilarityScore() {
		return score;
	}

	public String getScoreAlphaFile() {
		return scoreAlphaFile;
	}

	public void setScoreAlphaFile(String scoreAlphaFile) {
		this.scoreAlphaFile = scoreAlphaFile;
	}

	public String getExpectedProductFile() {
		return expectedProductFile;
	}

	public void setExpectedProductFile(String expectedProductFile) {
		this.expectedProductFile = expectedProductFile;
	}
}