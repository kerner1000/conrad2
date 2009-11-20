package calhoun.analysis.crf.executables.viewer;

import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.IOException;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.JScrollPane;
import javax.swing.JTable;

import calhoun.analysis.crf.Conrad;
import calhoun.analysis.crf.io.TrainingSequence;

public class ViterbiViewerApp extends JFrame {
	private static final long serialVersionUID = -5564046976928448201L;
	private static final int WINDOW_WIDTH = 800;
	private static final int WINDOW_HEIGHT = 600;
	ViterbiViewer recordViewer;

	public ViterbiViewerApp(String model, String data) throws IOException {
		super("Viterbi Path Viewer");
		addWindowListener(new WindowAdapter() {
			@Override
			public void windowClosing(WindowEvent e) {
				System.exit(0);
			}
		});
		
		Conrad r = Conrad.read(model);
		List<? extends TrainingSequence<?>> trainingData = r.getInputHandler().readTrainingData(data);
		
		recordViewer = new ViterbiViewer(r, trainingData.get(0));
		JScrollPane scrollPane = new JScrollPane(recordViewer);
		scrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
		setContentPane(scrollPane);
		//initializeMenuBar();
		setSize(WINDOW_WIDTH, WINDOW_HEIGHT);
	}

	public static void main(String[] args) throws Exception {
		String model = "";
		String data = "";
		if(args.length == 2){
			model = args[0];
			data = args[1];
		}
		final ViterbiViewerApp app = new ViterbiViewerApp(model, data);
		app.pack();
		app.setVisible(true);
	}

	static class ViterbiViewer extends JTable {
		private static final long serialVersionUID = 7592581578644968579L;
		ViterbiTableModel model; 

		public ViterbiViewer(Conrad crfModel, TrainingSequence seq) throws IOException {
			model = new ViterbiTableModel(crfModel, seq);
	        setModel(model);
	        setDefaultRenderer(String.class, model.new ViterbiCellRenderer());
	        //setShowGrid(true);
	        //setAutoCreateColumnsFromModel(true);
	        setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		}
	}
}
