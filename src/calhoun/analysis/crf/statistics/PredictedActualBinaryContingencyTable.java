package calhoun.analysis.crf.statistics;

import java.io.Serializable;
import java.text.DecimalFormat;
import java.text.NumberFormat;

import calhoun.util.Assert;

public class PredictedActualBinaryContingencyTable implements Serializable {
	private static final long serialVersionUID = 557844841355570852L;
	private int tp,fp,fn,tn;
	private boolean frozen;
	
	///////////////////////////////////////////
	// Constructors and modifiers are below:
	
	public PredictedActualBinaryContingencyTable() {
		this.tp = 0;
		this.fp = 0;
		this.fn = 0;
		this.tn = 0;
		frozen = false;
	}	

	public void set(int tp,int fp, int fn, int tn) {
		Assert.a(!frozen);
		this.tp = tp;
		this.fp = fp;
		this.fn = fn;
		this.tn = tn;
		frozen = true;
	}
	
	public void set(int tp,int fp, int fn) {
		Assert.a(!frozen);
		this.tp = tp;
		this.fp = fp;
		this.fn = fn;
		this.tn = -1;
		frozen = true;
	}
	
	public void forgetTN() {
		Assert.a(!frozen);
		Assert.a(tn ==0);
		tn = -1;
	}

	public void incrementTP() {
		Assert.a(!frozen);
		tp++;
	}
	
	public void incrementFP(){
		Assert.a(!frozen);
		fp++;
	}
	
	public void incrementFN(){
		Assert.a(!frozen);
		fn++;
	}
	
	public void incrementTN(){
		Assert.a(!frozen);
		Assert.a(tn>=0);
		tn++;
	}
	
	public void increment(boolean predicted, boolean actual) {
		Assert.a(!frozen);
		if ((predicted) && (actual)) { tp++; }
		if ((predicted) && (!actual)) { fp++; }
		if ((!predicted) && (actual)) { fn++; }
		if ((!predicted)&&(!actual)) {
			Assert.a(tn>=0);
			tn++;
		}
	}
	
	public void freeze() {
		frozen = true;
	}
		
	/////////////////////////////////////////////
	// Statistical measures
	
	public int getTP() {
		Assert.a(frozen);
		return tp;
	}
	
	public int getFP() {
		Assert.a(frozen);
		return fp;
	}
	
	public int getFN() {
		Assert.a(frozen);
		return fn;
	}
	
	public int getTN() {
		Assert.a(frozen);
		Assert.a(tn>=0);
		return tn;
	}
	
	/** Actual Positives */
	public int ap(){
		Assert.a(frozen);
		return (tp + fn);
	}

	/** Actual Negatives */
	public int an(){
		Assert.a(frozen);
		Assert.a(tn>=0);
		return (fp + tn);
	}
	
	/** Predicted Positives */
	public int pp(){
		Assert.a(frozen);
		return (tp + fp);
	}
	
	/** Predicted Negatives */
	public int pn(){
		Assert.a(frozen);
		Assert.a(tn>=0);
		return (tn+fn);
	}
	
	private boolean splitMargins() {
		if(tn>=0) { 
			if (tp+fn<=0) { return false; }
			if (tn+fp<=0) { return false; }
			if (tp+fp<=0) { return false; }
			if (tn+fn<=0) { return false; }
		} else {
			if (tp+fn<=0) { return false; }			
			if (tp+fp<=0) { return false; }
		}
		return true;
	}
	
	/** This is the Pearson correlation of two 0-1 random variables X=prediction and Y=reality 
	 *  CC = Cov(X,Y)/(Stddev(X)*Stddev(Y))
	 *  If either RV has zero variance, the CC is underfined: assertion faliure
	 *  If this contingency table is not tracking TN, then CC is undefined: assertion failure
	 */
	public double correlationCoefficient() {	
		Assert.a(frozen);
		Assert.a(splitMargins());
		double num = (tp*tn)-(fn*fp);
		double den2 = (double) (tp+fn)*(tn+fp)*(tp+fp)*(tn+fn);
		double cc = num/Math.sqrt(den2);
		return cc;
	}
	
	
	/** This is the average conditional probability.  Only defined if TN is being tracked, and
	 *  if all four marginal values are positive.
	 */
	public double averageConditionalProbability(){
		Assert.a(frozen);
		Assert.a(splitMargins());
		double acp = 0;
		acp += (double) tp/(tp+fn);
		acp += (double) tp/(tp+fp);
		acp += (double) tn/(tn+fp);
		acp += (double) tn/(tn+fn);
		acp /= 4.0;
		return acp;
	}
	
	/** This is the approximate correlation.  Only defined if TN is being tracked, and
	 *  if all four marginal values are positive.  Equal to
	 *  2*(averageConditionalProbability - 0.5)
	 */	
	public double approximateCorrelation(){
		Assert.a(frozen);
		return 2*(averageConditionalProbability() - 0.5);
	}
	
	/** This is the average of sensitivity and specifity.  This is equal to the limit of ACP as TN->infinity
	 *  This is defined even if we are not tracking TN.
	 */
	public double averageSensitivitySpecificity(){
		Assert.a(frozen);
		return 0.5*( sensitivity() + specificity() );
	}
	
	/** Sensitivity is TP/(TP+FN) is the fraction of actual events that are predicted.
	 *  If TP+FN is zero then sensitivity is undefined, resulting in assertion failure
	 *  Sensitivity is defined even if not tracking TN.
	 */
	public double sensitivity() {
		Assert.a(frozen);
		Assert.a(tp+fn > 0);
		return (double) tp/(tp+fn);
	}
	
	/** Specificity is the TP/(TP+FP) is the fraction of predicted events that are real.
	 *  If TP+FP=0 then specificity is undefined and result in assertion failure
	 *  Specificity is defined even if not tracking TN  
	 */
	public double specificity() {
		Assert.a(frozen);
		Assert.a(tp+fp > 0);
		return (double) tp/(tp+fp);
	}
	
	
	////////////////////////////////////////////////////
	// Printed summaries
	
	public String summarize() {
		String ret = "";
		NumberFormat d3 = new DecimalFormat();
		d3.setMinimumFractionDigits(3);
		d3.setMaximumFractionDigits(3);
		
		if (tn>=0) {
			ret += "( TP=" + tp + ", FP=" + fp + ", FN=" + fn + ", TN=" + tn + " ) ";
			ret += "( AP=" + ap() + ", AN=" + an() + ", PP=" + pp() + ", PN=" + pn() + " ) ";
			if (splitMargins()) {
				ret += "( CC=" + d3.format(correlationCoefficient()) + " ) ";
				ret += "( ACP=" + d3.format(averageConditionalProbability()) + " ) ";
				ret += "( AC=" + d3.format(approximateCorrelation()) + " ) ";
				ret += "( sens=" + d3.format(sensitivity()) + ", spec=" + d3.format(specificity()) + ", avSS=" + d3.format(averageSensitivitySpecificity()) + " ) ";
			} else {
				ret += "( margins not split ) ";
			}
		} else {
			ret += "( TP=" + tp + ", FP=" + fp + ", FN=" + fn + " ) ";			
			ret += "( AP=" + ap() + ", PP=" + pp() + " )";
			if (splitMargins() ) {
				ret += "( sens=" + d3.format(sensitivity()) + ", spec=" + d3.format(specificity()) + ", avSS=" + d3.format(averageSensitivitySpecificity()) + " ) ";
			} else {
				ret += "( margins not split ) ";
			} 
		}
		
		return ret;
	}

	
}
