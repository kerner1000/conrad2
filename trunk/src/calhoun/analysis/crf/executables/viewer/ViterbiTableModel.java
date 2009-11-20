package calhoun.analysis.crf.executables.viewer;

import java.awt.Color;
import java.awt.Component;

import javax.swing.JTable;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.DefaultTableCellRenderer;

import calhoun.analysis.crf.Conrad;
import calhoun.analysis.crf.CRFInference.InferenceResult;
import calhoun.analysis.crf.io.TrainingSequence;

public class ViterbiTableModel extends AbstractTableModel {
	private static final long serialVersionUID = 4315252264878822897L;
	Conrad crfModel;
	TrainingSequence seq;
	int[] viterbiPath;
	double[] bestScore;
	int nStates;

	public ViterbiTableModel(Conrad crfModel, TrainingSequence seq) {
		this.crfModel = crfModel;
		this.seq = seq;
		nStates = crfModel.getModel().getNumStates();
		InferenceResult result = crfModel.predict(seq);
		viterbiPath = result.hiddenStates;
		bestScore = null;
	}

	@Override
	public String getColumnName(int col) {
		return (col == 0) ? "Position" : "    " + Integer.toString(col) + "    ";
	}

	public int getRowCount() {
		return crfModel.getModel().getNumStates() + 2;
	}

	public int getColumnCount() {
		return seq.length() + 1;
	}

	@Override
	public Class getColumnClass(int col) {
		return String.class;
	}

	public Object getValueAt(int row, int col) {
		switch (row) {
		case 0:
			return col == 0 ? "Input" : seq.getInputSequence().getX(col - 1).toString();
		case 1:
			return col == 0 ? "Label" : Integer.toString(seq.getY(col - 1));
		default:
			return col == 0 ? crfModel.getModel().getStateName(row - 2) : Double.toString(bestScore[(row - 2) * nStates + col - 1]);
		}
	}

	@Override
	public boolean isCellEditable(int row, int col) {
		return false;
	}

	@Override
	public void setValueAt(Object value, int row, int col) {
		// this 1 line was screwing up the sort for almost a year.
		//m_tableSorter.fireTableDataChanged();
	}

	public class ViterbiCellRenderer extends DefaultTableCellRenderer {
		private static final long serialVersionUID = -1156070614840797662L;

		@Override
		public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
			super.getTableCellRendererComponent(table, value , isSelected, hasFocus, row, column);
			if(row > 1 && column > 0 && (row-2) == viterbiPath[column-1]) {
				setBackground(Color.RED);
			}
			else {
				setBackground(Color.WHITE);
			}
			return this;
		}
	}
}
