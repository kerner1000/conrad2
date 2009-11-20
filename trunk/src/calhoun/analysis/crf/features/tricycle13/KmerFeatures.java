package calhoun.analysis.crf.features.tricycle13;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.analysis.crf.AbstractFeatureManager;
import calhoun.analysis.crf.FeatureList;
import calhoun.analysis.crf.FeatureManagerNode;
import calhoun.analysis.crf.ModelManager;
import calhoun.analysis.crf.io.InputSequence;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.util.Assert;

/** trains on the data and then evaluates to P(state | label) for given Kmers.  Used for historical reasons.  Emission markov generally does better.  
 */
public class KmerFeatures extends AbstractFeatureManager<Character> implements FeatureManagerNode<Character> {
	private static final long serialVersionUID = 5959560033335736926L;
	private static final Log log = LogFactory.getLog(KmerFeatures.class);

	public enum Cardinality {
		SINGLE, PER_STATE, PER_KMER
	};

	static final int DEFAULT_RARE_THRESHOLD = 25;
	int startIx;
	ModelManager model;
	int nStates;
	List<int[]> kmerDefs;
	Map<String, double[]>[] counts;
	Map<String, Integer>[] kmerIds;
	Cardinality cardinality = Cardinality.PER_STATE;
	int rareThreshold = DEFAULT_RARE_THRESHOLD;

	public KmerFeatures(Cardinality cardinality) {
		this.cardinality = cardinality;
		this.kmerDefs = Collections.singletonList(new int[] { 0 });
	}

	public KmerFeatures() {
		this.kmerDefs = Collections.singletonList(new int[] { 0 });
	}

	public KmerFeatures(List<int[]> kmerDefs) {
		this.kmerDefs = kmerDefs;
	}

	public KmerFeatures(List<int[]> kmerDefs, Cardinality cardinality) {
		this.cardinality = cardinality;
		this.kmerDefs = kmerDefs;
	}

	public void setKmerDefinitions(List<List<Integer>> defs) {
		kmerDefs = new ArrayList<int[]>();
		for(List<Integer> def : defs) {
			int[] kmer = new int[def.size()];
			kmerDefs.add(kmer);
			for(int i=0; i<def.size(); ++i) {
				kmer[i] = def.get(i);
			}
		}
	}
	
	public void setRareThreshold(int threshold) {
		Assert.a(model == null, "Can't set threshold after training.");
		rareThreshold = threshold;
	}

	public int getNumFeatures() {
		switch (cardinality) {
		case SINGLE:
			return kmerDefs.size();
		case PER_STATE:
			return kmerDefs.size() * model.getNumStates();
		default:
			Assert.a(cardinality == Cardinality.PER_KMER);
			int num = 0;
			for (Map m : counts) {
				num += m.size();
			}
			return num;
		}
	}

	public String getFeatureName(int featureIndex) {
		int raw = featureIndex - startIx;
		int label = raw / kmerDefs.size();
		int kmer = raw % kmerDefs.size();
		String name = kmerName(kmer);
		String val = "";
		if (cardinality == Cardinality.PER_STATE) {
			val = model.getStateName(label) + ".";
		} else if (cardinality == Cardinality.PER_KMER) {
			val = "(InsertKMer)";
		}
		return "Kmer." + val + name;
	}

	transient InputSequence<? extends Character> lastSeq;
	transient int lastPos;
	String[] kmers;
	double[][] vals;

	public void evaluateNode(InputSequence<? extends Character> seq, int pos, int state, FeatureList result) {
		if (seq != lastSeq || pos != lastPos) {
			lastSeq = seq;
			lastPos = pos;
			for (int j = 0; j < kmerDefs.size(); ++j) {
				String kmer = getKmer(seq, pos, kmerDefs.get(j));
				vals[j] = counts[j].get(kmer);
				if (cardinality == Cardinality.PER_KMER) {
					kmers[j] = kmer;
				}
			}
		}
		for (int j = 0; j < kmerDefs.size(); ++j) {
			double[] kmerVals = vals[j];
			if (kmerVals == null || kmerVals[state] == 0.0)
				continue;
			int index;
			if (cardinality == Cardinality.SINGLE) {
				index = startIx + j;
			} else if (cardinality == Cardinality.PER_STATE) {
				index = startIx + j + kmerDefs.size() * state;
			} else {
				Assert.a(cardinality == Cardinality.PER_KMER);
				Integer id = kmerIds[j].get(kmers[j]);
				if (id == null)
					continue;
				index = id;
			}
			result.addFeature(index, kmerVals[state]);
		}
	}

	public String getKmer(InputSequence<? extends Character> seq, int pos, int[] def) {
		StringBuffer buf = new StringBuffer(def.length);
		for (int i = 0; i < def.length; ++i) {
			int loc = pos + def[i];
			if (loc < 0 || loc >= seq.length())
				return null;
			buf.append((Character) seq.getX(loc));
		}
		return buf.toString();
	}

	/** Returns a string representation of a given kmer definition */
	public String kmerName(int index) {
		int[] kmer = kmerDefs.get(index);
		StringBuffer ret = new StringBuffer();
		for (int i = 0; i < kmer.length; ++i) {
			if (i != 0) {
				ret.append(".");
			}
			ret.append(kmer[i]);
		}
		return ret.toString();
	}

	/** Returns an individual entry from the counts list. */
	public double getKmerProb(int kmerIndex, String kmer, int label) {
		return counts[kmerIndex].get(kmer)[label];
	}

	/** Computes the P(label | kmer) for each kmer across all of the training data. These will used as features values. */
	public void train(int startingIndex, ModelManager modelInfo, List<? extends TrainingSequence<? extends Character>> data) {
		startIx = startingIndex;
		model = modelInfo;
		nStates = model.getNumStates();
		vals = new double[kmerDefs.size()][];
		// Initialize an array of hash maps for each kmer
		counts = new HashMap[kmerDefs.size()];
		for (int i = 0; i < counts.length; ++i) {
			counts[i] = new HashMap();
		}
		if (cardinality == Cardinality.PER_KMER) {
			kmerIds = new HashMap[kmerDefs.size()];
			for (int i = 0; i < kmerIds.length; ++i) {
				kmerIds[i] = new HashMap();
			}
		}
		// Count the occurances of each kmer in each state
		for (TrainingSequence<? extends Character> seq : data) {
			int len = seq.length();
			for (int pos = 0; pos < len; ++pos) {
				int nKmers = kmerDefs.size();
				for (int i = 0; i < nKmers; ++i) {
					String kmer = getKmer(seq, pos, kmerDefs.get(i));
					if (kmer != null) {
						double[] val = (double[]) counts[i].get(kmer);
						if (val == null) {
							val = new double[nStates];
							Arrays.fill(val, 1);
							counts[i].put(kmer, val);
						}
						int state = seq.getY(pos);
						val[state] += 1;
					}
				}
			}
		}
		int kmerId = startIx;
		// Now compute probabilities for all kmers with at least KMER_THRESHOLD appearances
		for (int j = 0; j < kmerDefs.size(); ++j) {
			int kmerCounts = 0;
			for (Map.Entry<String, double[]> val : counts[j].entrySet()) {
				double[] cnts = val.getValue();
				double total = 0.0f;
				for (int i = 0; i < cnts.length; ++i) {
					total += cnts[i];
				}
				if (total >= rareThreshold) {
					kmerCounts += 1;
					if (cardinality == Cardinality.PER_KMER) {
						kmerIds[j].put(val.getKey(), kmerId);
						kmerId++;
					}
				}
				for (int i = 0; i < cnts.length; ++i) {
					// Just convert to P(label | kmer)
					if (total >= rareThreshold) {
						cnts[i] = Math.log(cnts[i] / total);
						Assert.a(!Double.isNaN(cnts[i]));
						Assert.a(!Double.isInfinite(cnts[i]));
					} else {
						cnts[i] = 0;
					}
				}
			}
			log.info(kmerCounts + " kmers of " + kmerName(j));
		}
	}
}
