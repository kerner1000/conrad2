package calhoun.analysis.crf.features.supporting;

import java.io.Serializable;
import java.util.List;

import calhoun.analysis.crf.io.InputSequence;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.seq.KmerHasher;
import calhoun.util.Assert;

public class MarkovPredictorLogprob implements Serializable {
	private static final long serialVersionUID = -3801967331270040642L;

	// For each state S, this class is capable of returning
	// log(p_S(x_i|x_i-h,...,x_i-2,x_i-1))    h = degree of predictor
	
	// Trained using a training sequence with hidden states labeled.
	// For training, one only wants to look at positions for which not only
	// the current hidden state but also the preceding h are as specified.
	
	// For example, for a fifth order predictor for the exon3 state, I might only
	// want to train using positions in the training data for which the preceding
	// five hidden states were (e1,e2,e3,12,e2).
	
	boolean trainedYet = false;
	List<int[]> history;
	int nStates;    // number of hidden states
	int[] kmerLengths;  // order of Markov model (5th order means x_i predicted based on previous five).
	KmerHasher[] hashers;  // One for each state (they might need different history lengths, so use different hashers)
	float[][][] logProb;  // Log probabilities.  First dimension is the state, second is indexed by the hash of the history, third is the hash of the current position.
	int maxLength;
	KmerHasher h; // for a single letter.  Used to index the third dimension of the log prob.
	
	int[] currentHash;
	char[] scratch;

	/** History is a list that contains one entry for each state.  Each entry in the list is in turn a list of the preceding states
	 * that we want to examine.
	 * 
	 * @param history
	 */
	public MarkovPredictorLogprob(List<int[]> history) {
		this.history = history;

		nStates = history.size();

		kmerLengths = new int[nStates];
		hashers = new KmerHasher[nStates];
		logProb = new float[nStates][][];
		currentHash = new int[nStates];

		h = new KmerHasher(KmerHasher.ACGTN, 1);
			
		maxLength = 0;
		for (int j=0; j<nStates; j++) {
			kmerLengths[j] = history.get(j).length;
			maxLength = Math.max(kmerLengths[j], maxLength);
			hashers[j] = new KmerHasher(KmerHasher.ACGTN, kmerLengths[j]);
			int historyRange = hashers[j].range();
			int range = h.range();
			logProb[j] = new float[historyRange][h.range()]; 
			for (int k=0; k<historyRange; k++) {
				for (int l=0; l<range; l++) {
					logProb[j][k][l] = 1.0f;   // Initializing the pseudocounts
				}
			}
		}
		scratch = new char[maxLength];
	}
	
	transient InputSequence<? extends Character> lastSeq = null;
	int lastPos = -1;
	char lastChar;
	int singleHash;
	public float logprob(int state, InputSequence<? extends Character> seq, int pos) {
		/* Returns the log probability of nucleotide at position pos in ISC given
		 * that the hidden state at position pos is state and given the previous several
		 * (about 4, depends on how trained) nucleotides in ISC.*/
		Assert.a(trainedYet);
		
		if (pos < maxLength || pos >= seq.length()) { 
			return (float) 0.0; 
		}		
		
		if(seq != lastSeq || pos != lastPos) {
			if(pos == lastPos + 1) {
				updateHashes(lastChar);
			}
			else {
				initHashes(seq, pos);
			}
			lastSeq = seq;
			lastPos = pos;
			lastChar = seq.getX(pos);
			singleHash = h.hash(lastChar);
		}
		
		return logProb[state][currentHash[state]][singleHash];
	}

	
	public void train(List<? extends TrainingSequence<? extends Character>> data) {
		if (trainedYet) { 
			return; 
		}
		
		// Loop through the data and increment all the counts
		for(TrainingSequence<? extends Character> seq : data) {
			int len = seq.length();

			initHashes(seq, maxLength);
			for (int pos=maxLength; pos<len; pos++) {
				int v = seq.getY(pos);
				boolean okHistory = true;
				for (int i=0; i<kmerLengths[v]; i++) {
					if (history.get(v)[i] != seq.getY(pos-kmerLengths[v]+i)) { 
						okHistory = false;
						break;
					}
				}
				
				char c = seq.getX(pos);
				if (okHistory) {
					logProb[v][currentHash[v]][h.hash(c)] += 1.0f;
				}
				updateHashes(c);
			}	
		}
		
		// Then normalize the counts
		for (int v=0; v<nStates; v++) {
			for (int j=0; j<logProb[v].length; j++) {
				float norm = (float) 0.0;
				for (int k=0; k<logProb[v][j].length; k++) {
					norm = norm + logProb[v][j][k];
				}
				double logNorm = Math.log(norm);
				for (int k=0; k<h.range(); k++) {
					logProb[v][j][k] = (float) (Math.log(logProb[v][j][k]) - logNorm);
				}				
			}
		}	
		trainedYet = true;
	}

	void initHashes(InputSequence<? extends Character> seq, int start) {
		// Fill the scratch buffer with the first entries, from 0 to maxLength-1
		for (int pos=0; pos<maxLength; pos++) {
			scratch[pos] = seq.getX(pos-maxLength+start);
			// Note: if you crash on above line with class cast exception, then a possible
			// cause is that you are using a composite input, but in you config file you
			// do not specify which input component this feature should be using.
		}
		for(int i = 0; i<currentHash.length; ++i) {
			currentHash[i] = hashers[i].hash(scratch, maxLength - kmerLengths[i]);
		}
	}

	void updateHashes(char c) {
		for(int i = 0; i<currentHash.length; ++i) {
			currentHash[i] = hashers[i].shiftHash(c, currentHash[i]);
		}
	}
}
