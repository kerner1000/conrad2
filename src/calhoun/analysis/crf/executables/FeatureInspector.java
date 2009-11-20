package calhoun.analysis.crf.executables;

import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.IOException;
import java.util.List;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextArea;
import javax.swing.JTextField;

import calhoun.analysis.crf.CRFObjectiveFunctionGradient;
import calhoun.analysis.crf.Conrad;
import calhoun.analysis.crf.io.TrainingSequence;
import calhoun.analysis.crf.solver.CacheProcessor;
import calhoun.analysis.crf.solver.CacheProcessorDeluxe;
import calhoun.analysis.crf.solver.MaximumLikelihoodGradient;
import calhoun.analysis.crf.solver.MaximumLikelihoodSemiMarkovGradient;
import calhoun.analysis.crf.solver.StandardOptimizer;
import calhoun.analysis.crf.solver.CacheProcessor.FeatureEvaluation;
import calhoun.analysis.crf.solver.CacheProcessor.LengthFeatureEvaluation;
import calhoun.analysis.crf.solver.CacheProcessor.SolverSetup;

/** provides a simple way to interactively examine feature evaluations.  
 * Takes a trained model and a data set on the command line and creates a window
 * where you can type in parameters to the evaluate function and view the results.
 */
public class FeatureInspector extends JFrame implements ActionListener {
	private static final long serialVersionUID = 4470742291972538174L;

	private static final int WINDOW_WIDTH = 800;
	private static final int WINDOW_HEIGHT = 600;

	CacheProcessorDeluxe cp;
	JTextField seq = new JTextField(5);
	JTextField start = new JTextField(5);
	JTextField end = new JTextField(5);
	JTextField from = new JTextField(5);
	JTextField to = new JTextField(5);
	JTextArea results = new JTextArea("Results:", 20, 40);
	
	public FeatureInspector() throws IOException {
		super("Feature Inspector");
		addWindowListener(new WindowAdapter() {
			@Override
			public void windowClosing(WindowEvent e) {
				System.exit(0);
			}
		});
		
		JPanel overall = new JPanel();
		overall.setLayout(new BorderLayout());
		JPanel pane = new JPanel();
		pane.add(new JLabel("Seq:"));
		pane.add(seq);
		pane.add(new JLabel("Start pos:"));
		pane.add(start);
		pane.add(new JLabel("End pos:"));
		pane.add(end);
		pane.add(new JLabel("From:"));
		pane.add(from);
		pane.add(new JLabel("To:"));
		pane.add(to);
		JButton button = new JButton("GetFeatures");
		button.addActionListener(this);
		pane.add(button);
		overall.add(pane, BorderLayout.NORTH);
		overall.add(results, BorderLayout.SOUTH);
		setContentPane(overall);
		//initializeMenuBar();
		setSize(WINDOW_WIDTH, WINDOW_HEIGHT);
	}

	public static void main(String[] args) throws Exception {
		if(args.length == 0) {
			args = new String[] { "test/working/interval13BaselineModelTest2.ser", "test/input/interval13/data/oneGeneTrain.interval13.txt"};
		}
		
		// First argument is the model
		Conrad c = Conrad.read(args[0]);
		
		// Second argument is the data
		List<? extends TrainingSequence<?>> data = c.getInputHandler().readTrainingData(args[1]);
		CRFObjectiveFunctionGradient crf = ((StandardOptimizer)c.getOptimizer()).getObjectiveFunction();
		CacheProcessor cp;
		if(crf instanceof MaximumLikelihoodSemiMarkovGradient)
			cp = ((MaximumLikelihoodSemiMarkovGradient) crf).getCacheProcessor();
		else
			cp = ((MaximumLikelihoodGradient) crf).getCacheProcessor();

		cp.setTrainingData(c.getModel(), data);
		
		final FeatureInspector app = new FeatureInspector();
		app.cp = (CacheProcessorDeluxe) cp;
		app.pack();
		app.setVisible(true);
	}

	public void actionPerformed(ActionEvent arg0) {
		results.setText(getMessage());
		invalidate();
	}
	
	String getMessage() {
		int seqNum = Integer.parseInt(seq.getText());
		int e = Integer.parseInt(end.getText());
		int t = Integer.parseInt(to.getText());
		int s = Integer.MIN_VALUE;
		if(start.getText().length()> 0)
			s = Integer.parseInt(start.getText());
		int f = Integer.MIN_VALUE;
		if(from.getText().length()> 0)
			f = Integer.parseInt(from.getText());

		SolverSetup solverSetup = cp.getSolverSetup();
		StringBuffer b = new StringBuffer();
		FeatureEvaluation eval;
		if(s == Integer.MIN_VALUE) {
			cp.evaluatePosition(seqNum, e);
			int potential;
			if(f == Integer.MIN_VALUE) {
				// Node Potential
				potential = t;
			}
			else {
				potential = solverSetup.transitionIndex.getQuick(f, t) + solverSetup.nStates;
				if(potential == -1 + solverSetup.nStates) {
					return "Transition is not in model.";
				}
			}
			eval = cp.getFeatureEvaluations()[potential];
			if(cp.invalidTransitions[(solverSetup.seqOffsets[seqNum]+e)*solverSetup.nPotentials + t])
				return "Invalid node";
			if(cp.invalidTransitions[(solverSetup.seqOffsets[seqNum]+e)*solverSetup.nPotentials + potential])
				return "Invalid edge";
			
		}
		else {
			cp.evaluateSegmentsEndingAt(seqNum, e);
			LengthFeatureEvaluation[][] l = cp.getLengthFeatureEvaluations();
			// Find the right state
			int stateIx = 0;
			while(solverSetup.statesWithLookback[stateIx].state != t) {
				++stateIx;
				if(stateIx == solverSetup.statesWithLookback.length)
					return "State "+t+" is not a semi-Markov state.";
			}
			// Find the right lookback
			int lookbackIx = 0;
			while(l[stateIx][lookbackIx].lookback != e-s) {
				if(l[stateIx][lookbackIx].lookback == -1)
					return "Lookback is not valid.";
				++lookbackIx;
			}
			b.append("Node features only:\n");
			eval = l[stateIx][lookbackIx].nodeEval;
		}

		int i = 0;
		while(eval.index[i] != -1) {
			b.append(eval.index[i]).append(": ").append(eval.value[i]).append("\n");
			++i;
		}
		if(i == 0)
			return "No features";
		return b.toString();
	}
}
