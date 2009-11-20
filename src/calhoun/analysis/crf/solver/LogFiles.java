package calhoun.analysis.crf.solver;

import java.io.BufferedWriter;

import calhoun.util.FileUtil;

public final class LogFiles {
	public String alphaFile = null;
	public String alphaLengthFile = null;
	public String betaLengthFile = null;
	public String expectFile = null;
	public String expectLengthFile = null;
	public String nodeMarginalFile = null;
	public String scoreAlphaFile = null;
	public String expectedProductFile = null;
	public String marginalsFile = null;
	public BufferedWriter alphaWriter = null;
	public BufferedWriter alphaLengthWriter = null;
	public BufferedWriter betaLengthWriter = null;
	public BufferedWriter expectWriter = null;
	public BufferedWriter expectLengthWriter = null;
	public BufferedWriter nodeMarginalWriter = null;
	public BufferedWriter scoreAlphaWriter = null;
	public BufferedWriter expectedProductWriter = null;
	public BufferedWriter marginalsWriter = null;
	
	public final void open() {
		alphaWriter = FileUtil.safeOpen(alphaFile);
		alphaLengthWriter = FileUtil.safeOpen(alphaLengthFile);
		betaLengthWriter = FileUtil.safeOpen(betaLengthFile);
		expectWriter = FileUtil.safeOpen(expectFile);
		expectLengthWriter = FileUtil.safeOpen(expectLengthFile);
		nodeMarginalWriter = FileUtil.safeOpen(nodeMarginalFile);
		scoreAlphaWriter = FileUtil.safeOpen(scoreAlphaFile);
		expectedProductWriter = FileUtil.safeOpen(expectedProductFile);
		marginalsWriter = FileUtil.safeOpen(marginalsFile);
	}

	public final void close() {
		FileUtil.safeClose(alphaWriter);
		FileUtil.safeClose(alphaLengthWriter);
		FileUtil.safeClose(betaLengthWriter);
		FileUtil.safeClose(expectWriter);
		FileUtil.safeClose(expectLengthWriter);
		FileUtil.safeClose(nodeMarginalWriter);
		FileUtil.safeClose(scoreAlphaWriter);
		FileUtil.safeClose(expectedProductWriter);
		FileUtil.safeClose(marginalsWriter);
	}
}
