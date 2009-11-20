package calhoun.analysis.crf.executables;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.geom.Ellipse2D;
import java.io.File;

import javax.swing.JPanel;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.renderer.category.StatisticalLineAndShapeRenderer;
import org.jfree.data.category.CategoryDataset;
import org.jfree.data.statistics.DefaultStatisticalCategoryDataset;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RefineryUtilities;

import calhoun.util.FileUtil;

public class JFreeChartUI2 extends ApplicationFrame {

	private static final long serialVersionUID = 2121829480061443584L;

	public JFreeChartUI2(JPanel chartPanel, String title) {
		super(title);
		chartPanel.setPreferredSize(new Dimension(1000, 700));
		setContentPane(chartPanel);
	}

	public static void main(String[] args) throws Exception {
		load("nuc");
		load("exon");
		load("gene");
	}
	
	static void load(String val) throws Exception {
		JPanel panel = loadChart(val);
		JFreeChartUI2 demo = new JFreeChartUI2(panel, val);
		demo.pack();
		RefineryUtilities.centerFrameOnScreen(demo);
		demo.setVisible(true);
	}

	static JPanel loadChart(String type) throws Exception {
		DefaultStatisticalCategoryDataset dataset = (DefaultStatisticalCategoryDataset) FileUtil.readObject(type+".dat");
		JFreeChart chart = createChart(dataset);
		ChartUtilities.saveChartAsPNG(new File(type+".png"), chart, 1000, 700);
		return new ChartPanel(chart);
	}
	
	private static JFreeChart createChart(CategoryDataset dataset) {

		// create the chart...
		JFreeChart chart = ChartFactory.createLineChart(null, // chart
																						// title
				"Training set size", // domain axis label
				"Accuracy", // range axis label
				dataset, // data
				PlotOrientation.VERTICAL, // orientation
				true, // include legend
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
		rangeAxis.setUpperBound(1.0);
		rangeAxis.setAutoRange(true);

		// customise the renderer...
		StatisticalLineAndShapeRenderer renderer = new StatisticalLineAndShapeRenderer();
        plot.setRenderer(renderer);
	    renderer.setShapesVisible(true);
	    renderer.setDrawOutlines(true);
	    renderer.setUseFillPaint(true);
	    renderer.setFillPaint(Color.white);
//	    renderer.setSeriesVisible(false);
//	    renderer.setShapesVisible(false);
//	    renderer.setLinesVisible(false);
//	    renderer.setSeriesVisible(0, true);
//	    renderer.setSeriesVisible(1, true);
//	    renderer.setSeriesVisible(2, true);
//	    renderer.setSeriesVisible(3, true);
//	    renderer.setSeriesVisible(4, true);
//	    renderer.setSeriesVisible(5, false);
//	    renderer.setSeriesVisible(6, false);
//	    renderer.setSeriesVisible(7, false);
//	    renderer.setSeriesVisible(8, false);
//	    renderer.setSeriesVisible(9, false);
//	    renderer.setSeriesVisible(10, false);
//	    renderer.setSeriesVisible(11, false);
//	    renderer.setSeriesVisible(12, false);
	    renderer.setSeriesStroke(0, new BasicStroke(3.0f));
	    renderer.setSeriesOutlineStroke(0, new BasicStroke(2.0f));
	    renderer.setSeriesShape(0, new Ellipse2D.Double(-5.0, -5.0, 10.0, 10.0));
	    return chart;
	}
}