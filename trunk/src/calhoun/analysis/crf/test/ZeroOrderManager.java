package calhoun.analysis.crf.test;

import calhoun.analysis.crf.Conrad;
import calhoun.analysis.crf.features.generic.EndFeatures;
import calhoun.analysis.crf.features.generic.IndicatorEdges;
import calhoun.analysis.crf.features.generic.StartFeatures;
import calhoun.analysis.crf.features.tricycle13.KmerFeatures;
import calhoun.analysis.crf.io.InputHandlerInterleaved;
import calhoun.analysis.crf.io.OutputHandlerGeneCallStats;
import calhoun.analysis.crf.io.StringInput;
import calhoun.analysis.crf.solver.MaximumLikelihoodGradient;
import calhoun.analysis.crf.solver.StandardOptimizer;
import calhoun.analysis.crf.solver.Viterbi;


public class ZeroOrderManager extends ZeroOrderModel {
	private static final long serialVersionUID = 3959312826759045449L;

	public ZeroOrderManager() {
		addFeatureManager(new StartFeatures());
		addFeatureManager(new EndFeatures());
		addFeatureManager(new KmerFeatures());
		addFeatureManager(new IndicatorEdges());
	}

	public static Conrad getCRF() {
		Conrad ret = new Conrad();
		ret.setInference(new Viterbi());
		ret.setModel(new ZeroOrderManager());
		ret.setInputHandler(new InputHandlerInterleaved(new StringInput()));
		OutputHandlerGeneCallStats stats = new OutputHandlerGeneCallStats(ret.getModel(), ret.getInputHandler());
		stats.setWriteTrainingData(true);
		ret.setOutputHandler(stats);
		StandardOptimizer opt = new StandardOptimizer();
		opt.setDebugLevel(2);
		opt.setRequireConvergence(false);
		opt.setEpsForConvergence(0.0000005);
		opt.setObjectiveFunction(new MaximumLikelihoodGradient());
		ret.setOptimizer(opt);
		return ret;
	}
}
