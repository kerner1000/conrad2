package calhoun.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.Reader;
import java.io.StringWriter;
import java.util.ArrayList;

import bmsi.util.Diff;
import bmsi.util.DiffPrint;

/**
 * Computes the 'diff' difference between two files. By default, whitespace in the files is ignored.
 */
public class FileDiff
{
	private String[] firstLines;
	private String[] secondLines;

	private boolean closeReaders;
	private Reader firstLinesReader;
	private Reader secondLinesReader;
	private boolean ignoreWhitespace = true;
	private String commentLineIndicator;
	private Diff.change change = null;

	public static boolean filesMatch(String first, String second) throws IOException {
		return !new FileDiff(first, second).execute();
	}

	/**	Constructor takes 2 file names - opens the files using UTF-8 encoding.
	*/
	public FileDiff(String first, String second) throws IOException
	{
		this(new InputStreamReader( new FileInputStream(first), "UTF-8" ), new InputStreamReader( new FileInputStream(second),"UTF-8"), true );
	}

	/**	Constructor takes 2 Files - opens the files using UTF-8 encoding.
	*/
	public FileDiff(File first, File second) throws IOException
	{
		this(new InputStreamReader( new FileInputStream(first), "UTF-8" ), new InputStreamReader( new FileInputStream(second),"UTF-8"), true );
	}

	public FileDiff(Reader first, Reader second)
	{
		this( first,second, false );
	}

	public FileDiff(Reader first, Reader second, boolean closeReaders )
	{
		this.firstLinesReader = first;
		this.secondLinesReader = second;
		this.closeReaders = closeReaders;
	}

	/**
	 * Sets whether to ignore whitespace. By default, whitespace is ignored.
	 * @see #getIgnoreWhitespace
	 */
	public void setIgnoreWhitespace(boolean b)
	{
		ignoreWhitespace = b;
	}

	/**
	 * @see #setIgnoreWhitespace
	 */
	public boolean getIgnoreWhitespace()
	{
		return ignoreWhitespace;
	}

	/**
	 * Sets the comment line indicator string.  If set, lines beginning with this string are not used in the diff comparison.
	 * @see #getCommentLineIndicator
	 */
	public void setCommentLineIndicator(String commentLineIndicator)
	{
		this.commentLineIndicator = commentLineIndicator;
	}

	/**
	 * Gets the comment line indicator string.
	 * @see #setCommentLineIndicator
	 */
	public String getCommentLineIndicator()
	{
		return commentLineIndicator;
	}

	/**
	 * @return whether the two files were different.
	 * @see #toString
	 */
	public boolean execute() throws IOException
	{
		firstLines = readLines(firstLinesReader);
		secondLines = readLines(secondLinesReader);
		if ( closeReaders )
		{
			firstLinesReader.close();
			secondLinesReader.close();
		}
		Diff diff = new Diff(firstLines, secondLines);
		change = diff.diff_2(false);
		return change != null;
	}

	/**
	 * @return a diff-style listing of the differences between the 2 files. If the files were identical,
	 * returns the empty string ('').
	 */
	@Override
	public String toString()
	{
		if ( change != null )
		{
			DiffPrint.NormalPrint np = new DiffPrint.NormalPrint(firstLines, secondLines);
			StringWriter str = new StringWriter();
			PrintWriter writer = new PrintWriter(str);
			np.setWriter(writer);
			np.print_script(change);
			writer.close();
			return str.toString();
		}
		else
		{
			return "";
		}
	}

	private String[] readLines(Reader r) throws IOException
	{
		BufferedReader reader = new BufferedReader(r);
		ArrayList lines = new ArrayList();
		String line;
		while ( ( line = reader.readLine() ) != null )
		{
			if ( ignoreWhitespace )
				line = line.trim();

			if ( (!ignoreWhitespace || !"".equals(line)) &&
				  (commentLineIndicator == null || !line.startsWith(commentLineIndicator)) )
			{
				lines.add(line);
			}
		}
		String[] result = new String[lines.size()];
		return (String[])lines.toArray(result);
	}
}
