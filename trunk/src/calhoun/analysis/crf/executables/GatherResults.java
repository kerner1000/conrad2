package calhoun.analysis.crf.executables;

import org.apache.commons.math.stat.StatUtils;
import org.apache.commons.math.stat.descriptive.moment.StandardDeviation;
import org.jfree.data.statistics.DefaultStatisticalCategoryDataset;

import calhoun.analysis.crf.io.OutputHandlerGeneCallStats;
import calhoun.analysis.crf.statistics.PredictedActualBinaryContingencyTable;
import calhoun.util.ColtUtil;
import calhoun.util.FileUtil;

/**
 * Given a list of models, create several data sets to calculate: Average
 * statistics for accuracy with the 800 gene training sets: Nucleotide
 * sensitivity/specificity Exon sens./spec. Gene sens./spec.
 * 
 * Nucleotide accuracy graph across training set sizes Gene accuracy across
 * training set sizes
 * 
 */
public class GatherResults {
	static String[] graphList = new String[] { "crypto_mean","crypto_2way_gen","crypto_2way_cml","crypto_2way_mea","crypto_2way_mean","crypto_2way_gaps_foot_mean","crypto_5way_mean","crypto_5way_gaps_foot_mean","crypto_5way_gaps_foot_halfest_mean" };
	//static String[] specList = new String[] { "CND1_800_50/crypto_2way_cml_CND1_800_50","CND1_800_50/crypto_2way_gen_CND1_800_50" };
	static String[] specList = new String[] { "crypto_mean","crypto_2way_gen","crypto_2way_cml","crypto_2way_mea","crypto_2way_mean","crypto_2way_gaps_foot_mean","crypto_5way_mean","crypto_5way_gaps_foot_mean","crypto_5way_gaps_foot_halfest_mean" };
	//static String[] specList = new String[] { "crypto_2way_gen", "crypto_2way_cml", "crypto_2way_mea", "crypto_2way_mean", "crypto_5way_cml", "crypto_5way_mea", "crypto_5way_mean", "crypto_5way_gaps_foot_cml", "crypto_5way_gaps_foot_mea", "crypto_5way_gaps_foot_mean" };
	static String[] speciesList = new String[] { "AN1_A_oryzae_RIB40", "AN1_A_oryzae_RIB40_AFL2", "AN1_AC2_NF2_AF2_ANIG1_AT1_AFL2_A_oryzae_RIB40","AN1_AC2_NF2_AF2_ANIG1_AT1_AFL2_A_oryzae_RIB40_CI2", "AN1_AF2", "AN1_AF2_A_oryzae_RIB40", "AN1_AF2_A_oryzae_RIB40_AFL2", "AN1_AF2_A_oryzae_RIB40_CI2", "AN1_AF2_ANIG1_A_oryzae", "AN1_CI2" };
	static String[] modelList = new String[] { "comp_gaps_foot_mean" };
/*	static String[] graphList = new String[] { "baseline_2way", "baseline_2way_untied_cml",
			"baseline_2way_untied_all_aof", "baseline_2way_untied_all_aof_nucsplice",
			"baseline_2way_untied_all_aof_exonbound", "baseline_untied_all_cml", "baseline_5way_untied_cml",
			"baseline_2way_gaps_foot_untied_all_cml", "baseline_5way_gaps_foot_untied_all_cml",
			"baseline_5way_gaps_foot_blast_untied_all_cml", "baseline_5way_gaps_foot_blast_untied_all_reg_cml", "baseline_5way_gaps_foot_blast_est_untied_all_cml" };*/
//	static String[] graphList = new String[] { "baseline_2way_untied_cml", "baseline_5way_untied_cml" };
//	static String[] graphList = new String[] { "baseline_5way_gaps_foot_blast_untied_all_reg_cml" };
	
	static String[] sizes = new String[] { "50", "100", "200", "400", "600", "800", "1000" };
	static String template="y:/scripts/conrad/daved/results/";
	static String specSet= "CND1_600";

	// baseline_2way_untied_all_aof_nucsplice,baseline_2way_untied_cml,baseline_2way_gaps_foot_untied_all_cml,baseline_5way_gaps_foot_blast_est_untied_all_cml,baseline_5way_gaps_foot_blast_untied_all_cml
	enum Stat {
		NUC_SENS(0), NUC_SPEC(1), EXON_SENS(2), EXON_SPEC(3), GENE_SENS(4), GENE_SPEC(5);
		Stat(int value) {
			this.value = value;
		}

		int value;
	};

	public static void main(String[] args) throws Exception {
		if (args[0].equals("spec")) {
			//System.out.print("\tNuc. sens.\tNuc. spec.\tExon sens.\tExon spec.\tGene sens.\tGene spec.\n");
			String ret = "";
			for (String model : specList) {
				ret += sensSpec(template, specSet+"/"+model+"_"+specSet);
			}
			FileUtil.writeFile("sensSpec.txt", ret);
		} else if (args[0].equals("graph")) {
			graphs(template, graphList);
		} else if (args[0].equals("asp")) {
			String ret="\tNuc. sens.\tNuc. spec.\tExon sens.\tExon spec.\tGene sens.\tGene spec.\n";
			ret += sensSpec(template, speciesList[0]+"/ss_mean_"+speciesList[0]);
			for (String species : speciesList) {
				for (String model : modelList) {
					ret += sensSpec(template, species+"/"+model+"_"+species);
				}
			}
			FileUtil.writeFile("aspSensSpec.txt", ret);
		}
	}

	public static void graphs(String dir, String[] graphList) throws Exception {
		DefaultStatisticalCategoryDataset nuc = new DefaultStatisticalCategoryDataset();
		DefaultStatisticalCategoryDataset exon = new DefaultStatisticalCategoryDataset();
		DefaultStatisticalCategoryDataset gene = new DefaultStatisticalCategoryDataset();

		// Go through each file in the graph list and generate every training
		// set size
		StandardDeviation s = new StandardDeviation();
		for (String baseModel : graphList) {
			for (String size : sizes) {
				double[] nucWeights = new double[10];
				double[] exonWeights = new double[10];
				double[] geneWeights = new double[10];
				double[] nucWeightsT = new double[10];
				double[] exonWeightsT = new double[10];
				double[] geneWeightsT = new double[10];
				for (int i = 0; i < 10; ++i) {
					String file = dir+"CND1_"+size+"/"+baseModel+"_CND1_"+size+"_"+i+".ser.test_"+i+".dat";
					try {
						OutputHandlerGeneCallStats.Results res = (OutputHandlerGeneCallStats.Results) FileUtil
								.readObject(file);
						nucWeights[i] = res.ctCodingNucleotide.sensitivity();
						exonWeights[i] = res.ctExons.sensitivity();
						float total = res.perfect + res.imperfect;
						geneWeights[i] = res.perfect / total;
						System.out.println(String.format("%s %.1f", file, geneWeights[i]*100));

						file = dir+"CND1_"+size+"/"+baseModel+"_CND1_"+size+"_"+i+".ser.train_"+i+".dat";
						res = (OutputHandlerGeneCallStats.Results) FileUtil.readObject(file);
						nucWeightsT[i] = res.ctCodingNucleotide.sensitivity();
						exonWeightsT[i] = res.ctExons.sensitivity();
						total = res.perfect + res.imperfect;
						geneWeightsT[i] = res.perfect / total;
						//System.out.println(geneWeightsT[i] + file);
					} catch (Exception ex) {
						System.out.println(String.format("Skipping %s", file));
						continue;
					}
				}
				nuc.add(myMean(nucWeights), s.evaluate(filter(nucWeights)), baseModel, size);
				exon.add(myMean(exonWeights), s.evaluate(filter(exonWeights)), baseModel, size);
				gene.add(myMean(geneWeights), s.evaluate(filter(geneWeights)), baseModel, size);
				nuc.add(myMean(nucWeightsT), s.evaluate(filter(nucWeightsT)), baseModel + "_training", size);
				exon.add(myMean(exonWeightsT), s.evaluate(filter(exonWeightsT)), baseModel + "_training", size);
				gene.add(myMean(geneWeightsT), s.evaluate(filter(geneWeightsT)), baseModel + "_training", size);
			}
		}
		// For each training set size, average the 10 replicates to create a
		// data point
		// Data points exist for Nucleotide accuracy, exon accuracy, and gene
		// accuracy (sensitivity)
		FileUtil.writeObject("nuc.dat", nuc);
		FileUtil.writeObject("exon.dat", exon);
		FileUtil.writeObject("gene.dat", gene);
		FileUtil.writeFile("nuc.txt", printSet(nuc));
		FileUtil.writeFile("exon.txt", printSet(exon));
		FileUtil.writeFile("gene.txt", printSet(gene));
	}

	static double myMean(double[] vals) {
		return StatUtils.mean(filter(vals));
	}

	static double[] filter(double[] vals) {
		int count = 0;
		for(double val : vals) {
			if(val != 0.0) {
				count++;
			}
		}
		double[] nonZero = new double[count];
		count = 0;
		for(double val : vals) {
			if(val != 0.0) {
				nonZero[count++] = val;
			}
		}
		return nonZero;
	}
	
	static String printSet(DefaultStatisticalCategoryDataset set) {
		StringBuffer buf = new StringBuffer();
		for(int col=0; col<set.getRowCount(); ++col) {
			buf.append("\t").append(set.getRowKey(col));
		}
		for(int col=0; col<set.getRowCount(); ++col) {
			buf.append("\t").append(set.getRowKey(col)+"StdDev");
		}
		buf.append("\n");
		for(int row = 0; row < set.getColumnCount(); ++row) {
			buf.append(set.getColumnKey(row));
			for(int col=0; col<set.getRowCount(); ++col) {
				buf.append("\t").append(set.getMeanValue(col, row).floatValue()*100);
			}
			for(int col=0; col<set.getRowCount(); ++col) {
				buf.append("\t").append(set.getStdDevValue(col, row).floatValue()*100);
			}
			buf.append("\n");
		}
		return buf.toString();
	}
	
	public static String sensSpec(String dir, String model) throws Exception {
		double[][] values = new double[Stat.values().length][10];
		// Iterate through the 10 replicates for the model to pull stats
		System.out.println(model);
		for (int i = 0; i < 10; ++i) {
			try {
				String file = dir+model+"_"+i+".ser.test_"+i+".dat";
				OutputHandlerGeneCallStats.Results res = (OutputHandlerGeneCallStats.Results) FileUtil.readObject(file);
				values[Stat.NUC_SENS.value][i] = res.ctCodingNucleotide.sensitivity();
				values[Stat.NUC_SPEC.value][i] = res.ctCodingNucleotide.specificity();
				values[Stat.EXON_SENS.value][i] = res.ctExons.sensitivity();
				values[Stat.EXON_SPEC.value][i] = res.ctExons.specificity();
				float total = res.perfect + res.imperfect;
				values[Stat.GENE_SENS.value][i] = res.perfect / total;
				for (int k = 0; k < 6; ++k) {
					PredictedActualBinaryContingencyTable tab = res.ctTransitions.get(1 + k);
					total += tab.getFP() - tab.getFN();
				}
				values[Stat.GENE_SPEC.value][i] = res.perfect / total;
			} catch (Exception ex) {
				System.out.println(String.format("Skipping %s %d", model, i));
				continue;
			}
		}

		// Average the stats
		StandardDeviation s = new StandardDeviation();
		DefaultStatisticalCategoryDataset dataset = new DefaultStatisticalCategoryDataset();
		StringBuffer b = new StringBuffer();

		b.append(model);
		for (int j = 0; j < 6; ++j) {
			System.out.println(Stat.values()[j]+" "+ColtUtil.format(values[j]));
			double mean = myMean(values[j]);
			double stddev = s.evaluate(filter(values[j]));
			dataset.add(mean, stddev, model, Stat.values()[j].name());
			//b.append("\t").append(String.format("%.2f", mean * 100)); // .append("\t").append(String.format("%.2f",
																		// stddev*100));
			b.append("\t").append(String.format("%.1f\t%.1f", mean * 100, stddev*100));
		}
		b.append("\n");
		// System.out.println(model);
		//FileUtil.writeFile("sensSpec.txt", b.toString());
		//FileUtil.writeObject("chart.dat", dataset);
		return b.toString();
	}
}
