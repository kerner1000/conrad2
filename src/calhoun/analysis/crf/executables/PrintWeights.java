package calhoun.analysis.crf.executables;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math.stat.StatUtils;
import org.apache.commons.math.stat.descriptive.moment.StandardDeviation;
import org.jfree.data.statistics.DefaultStatisticalCategoryDataset;

import calhoun.analysis.crf.Conrad;
import calhoun.util.Assert;
import calhoun.util.FileUtil;

public class PrintWeights {
	// This one shows that mean reduces standard dev over mea.  cml is similar but puts lower weight on exon lengths
	//static String[] models = new String[] { "CND1_600/crypto_2way_cml_CND1_600","CND1_600/crypto_2way_mea_CND1_600","CND1_600/crypto_2way_mean_CND1_600" };
	// Aspergillus and crypto show similar patterns
	//static String[] models = new String[] { "CND1_600/crypto_2way_mean_CND1_600", "AN1_AF2/ss_mean_AN1_AF2", "AN1_AF2_A_oryzae_RIB40/ss_mean_AN1_AF2_A_oryzae_RIB40" };
	// Adding in gaps and footprint is reasonable
	//static String[] models = new String[] { "CND1_600/crypto_2way_mean_CND1_600", "CND1_600/crypto_2way_gaps_foot_mean_CND1_600" };
	// Just a single model
	//static String[] models = new String[] { "CND1_800/crypto_2way_mean_CND1_800" };
	// Just a single model
	//static String[] models = new String[] { "CND1_600/crypto_2way_mean_CND1_600", "AN1_AF2/comp_gaps_foot_mean_AN1_AF2" };
	// For paper
	static String[] models = new String[] { "CND1_600/crypto_2way_mean_CND1_600", "AN1_AF2/comp_mean_AN1_AF2" };
	// EST is crap
	//static String[] models = new String[] { "CND1_600/crypto_5way_gaps_foot_blast_est_mean_CND1_600","CND1_600/crypto_5way_gaps_foot_blast_halfest_mean_CND1_600" };
	// half-EST6
	//static String[] models = new String[] { "CND1_600/crypto_2way_gaps_foot_mean_CND1_600","CND1_600/crypto_5way_gaps_foot_blast_halfest_mean_CND1_600" };
	//static String[] models = new String[] { "CND1_600/crypto_2way_cml_CND1_600","CND1_600/crypto_2way_mea_CND1_600","CND1_600/crypto_2way_mean_CND1_600","CND1_600/crypto_5way_gaps_foot_mean_CND1_600","CND1_600/crypto_5way_gaps_foot_blast_est_mean_CND1_600","CND1_600/crypto_5way_gaps_foot_blast_halfest_mean_CND1_600" };
	//static String[] models = new String[] { "CND1_400/crypto_2way_cml_CND1_400","CND1_400/crypto_2way_mea_CND1_400","CND1_400/crypto_2way_mean_CND1_400" };
	//static String[] models = new String[] { "AN1_AF2/ss_mean_AN1_AF2", "AN1_AF2/comp_gaps_foot_mean_AN1_AF2", "AN1_AF2_A_oryzae_RIB40/comp_gaps_foot_mean_AN1_AF2_A_oryzae_RIB40" };
	static String template="y:/scripts/conrad/daved/results/";
	
	public static void main(String[] args) throws Exception {
		DefaultStatisticalCategoryDataset dataset = new DefaultStatisticalCategoryDataset();
		for(String model : models) {
			chart(dataset, template, model);
		}
		FileUtil.writeObject("Chart.dat", dataset);
		FileUtil.writeFile("weights.txt", printSet(dataset));
	}
	
	static String printSet(DefaultStatisticalCategoryDataset set) {
		StringBuffer buf = new StringBuffer();
		for(int col=0; col<set.getRowCount(); ++col) {
			if(col != 0)
				buf.append("\t");
			buf.append(set.getRowKey(col));
		}
		for(int col=0; col<set.getRowCount(); ++col) {
			buf.append("\t").append(set.getRowKey(col)+"StdDev");
		}
		buf.append("\n");
		for(int row = 0; row < set.getColumnCount(); ++row) {
			buf.append(set.getColumnKey(row));
			for(int col=0; col<set.getRowCount(); ++col) {
				buf.append("\t").append(safePrint(set.getMeanValue(col, row)));
			}
			for(int col=0; col<set.getRowCount(); ++col) {
				buf.append("\t").append(safePrint(set.getStdDevValue(col, row)));
			}
			buf.append("\n");
		}
		return buf.toString();
	}
	
	static String safePrint(Object o) {
		return (o == null) ? "" : o.toString();
	}
	
	public static void chart(DefaultStatisticalCategoryDataset dataset, String template, String model) throws Exception {
		int REP=10;
		double[][] weights = new double[REP][];

		Conrad r = null;
		int len = -1;
		for (int i = 0; i < REP; ++i) {
			String file = template+model+"_" + i + ".ser";
			System.out.println(file);
			try {
				r = Conrad.read(file);
				weights[i] = r.getWeights();
				len = weights[i].length;
			}
			catch(IOException ex) {
				ex.printStackTrace();
				System.out.println("Skipping "+file+" "+ex.getMessage());
				weights[i] = new double[len];
			}
		}
		List<String> names = new ArrayList();
		for (int i = 0; i < r.getNumFeatures(); ++i) {
			names.add(r.getFeatureName(i));
		}

		DefaultStatisticalCategoryDataset localDataset = new DefaultStatisticalCategoryDataset();

		Assert.a(r.getNumFeatures() == weights[REP-1].length, r.getNumFeatures()," ", weights[REP-1].length);
		Assert.a(names.size() == weights[REP-1].length, names.size()," ", weights[REP-1].length);
		StringBuffer b = new StringBuffer();
		double[] featureWeights = new double[REP];
		System.out.println(names.size());
		System.out.println(weights[REP-1].length);
		for (int j = 0; j < names.size(); ++j) {
			b.append(names.get(j));
			for (int i = 0; i < REP; ++i) {
				featureWeights[i] = weights[i][j];
			}
			b.append("\t").append(StatUtils.mean(featureWeights));
			StandardDeviation s = new StandardDeviation();
			b.append("\t").append(s.evaluate(featureWeights));
			dataset.add(StatUtils.mean(featureWeights), s.evaluate(featureWeights), model, names.get(j));
			localDataset.add(StatUtils.mean(featureWeights), s.evaluate(featureWeights), model, names.get(j));

			for (int i = 0; i < REP; ++i) {
				b.append("\t").append(featureWeights[i]);
			}
			b.append("\n");
		}
		System.out.print(b.toString());
		//FileUtil.writeObject(model+"Chart.dat", localDataset);
	}
}
