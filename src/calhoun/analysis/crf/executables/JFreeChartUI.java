package calhoun.analysis.crf.executables;

import java.awt.Color;
import java.awt.Dimension;
import java.io.File;
import java.text.DecimalFormat;
import java.text.FieldPosition;

import javax.swing.JPanel;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.labels.ItemLabelAnchor;
import org.jfree.chart.labels.ItemLabelPosition;
import org.jfree.chart.labels.StandardCategoryItemLabelGenerator;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.renderer.category.StatisticalBarRenderer;
import org.jfree.data.category.CategoryDataset;
import org.jfree.data.statistics.DefaultStatisticalCategoryDataset;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RefineryUtilities;
import org.jfree.ui.TextAnchor;

import calhoun.util.FileUtil;

public class JFreeChartUI extends ApplicationFrame {

	private static final long serialVersionUID = 3968640935099256212L;
	public JFreeChartUI(JPanel chartPanel, String title) {
		super(title);
		chartPanel.setPreferredSize(new Dimension(1000, 700));
		setContentPane(chartPanel);
	}

	//static String[] models = new String[] { "asp_gaps_foot_untied_all_cml", "baseline_2way_untied_cml", "baseline_5way_gaps_foot_blast_untied_all_cml", "baseline_5way_gaps_foot_blast_untied_all_reg_cml", "baseline_5way_gaps_foot_blast_est_untied_all_cml" };
	//static String[] models = new String[] { "baseline_5way_gaps_foot_blast_untied_all_reg_cml"};
	static String[] models = new String[] { ""};
	//static String[] models = new String[] { "baseline_untied_all_cml", "baseline_2way_untied_cml", "baseline_5way_untied_loweps_cml" };

	public static void main(String[] args) throws Exception {
//		JPanel panel = new JPanel();
//		for(String model : models) {
//			panel.add(createChart(model));
//		}
//		JFreeChartUI demo = new JFreeChartUI(panel, "Statistical Bar Chart Demo");
		for(String model : models) {
			JFreeChartUI demo = new JFreeChartUI(createChart(model), model);
			demo.pack();
			RefineryUtilities.centerFrameOnScreen(demo);
			demo.setVisible(true);
		}
	}

	static ChartPanel createChart(String model) throws Exception {
		DefaultStatisticalCategoryDataset dataset = (DefaultStatisticalCategoryDataset) FileUtil.readObject(model+"Chart.dat");
		JFreeChart chart = createChart(dataset);
		ChartUtilities.saveChartAsPNG(new File(model+"Weights.png"), chart, 1000, 700);
		return new ChartPanel(chart);
	}
	
	static class NumberFormatMine extends DecimalFormat {
		private static final long serialVersionUID = -2978723914826541314L;

		@Override
		public StringBuffer format(double arg0, StringBuffer arg1, FieldPosition arg2) {
			StringBuffer ret = super.format(arg0, arg1, arg2);
			System.out.println(ret);
			return ret;
		}
	}
	
	private static JFreeChart createChart(CategoryDataset dataset) {

		// create the chart...
		JFreeChart chart = ChartFactory.createLineChart(null, // chart
																						// title
				"Features", // domain axis label
				"Weights", // range axis label
				dataset, // data
				PlotOrientation.HORIZONTAL, // orientation
				false, // include legend
				true, // tooltips
				false // urls
				);

		chart.setBackgroundPaint(Color.white);
		chart.setBorderVisible(false);
		//chart.getCategoryPlot().
		
		CategoryPlot plot = (CategoryPlot) chart.getPlot();
		plot.setBackgroundPaint(Color.white);
		plot.setRangeGridlinePaint(Color.GRAY);
		//StandardCategoryItemLabelGenerator gen = new StandardCategoryItemLabelGenerator("%f", new NumberFormatMine());
		//plot.getRenderer().setBaseItemLabelGenerator(gen);
		
		// customise the range axis...
		NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
		rangeAxis.setAxisLineVisible(false);
		rangeAxis.setStandardTickUnits(NumberAxis.createStandardTickUnits());
		rangeAxis.setAutoRangeIncludesZero(false);
		rangeAxis.setNumberFormatOverride(new NumberFormatMine());

		// customise the renderer...
		StatisticalBarRenderer renderer = new StatisticalBarRenderer();
		renderer.setBaseFillPaint(Color.blue);
		renderer.setErrorIndicatorPaint(Color.black);
		renderer.setIncludeBaseInRange(false);
		plot.setRenderer(renderer);

		renderer.setItemLabelGenerator(new StandardCategoryItemLabelGenerator("{2}", new NumberFormatMine()));
		renderer.setBaseItemLabelsVisible(true);
		renderer.setPositiveItemLabelPosition(new ItemLabelPosition(ItemLabelAnchor.INSIDE6, TextAnchor.BOTTOM_CENTER));
		// OPTIONAL CUSTOMISATION COMPLETED.
		return chart;
	}
}