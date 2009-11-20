package calhoun.seq;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.util.ConfigException;
import calhoun.util.FileUtil;

/** Used for writing fasta files.  Formats sequences to a specific line length.
 * Allows multiple sequences to be to the same file.
 */
public class FastaWriter {
	private static final Log log = LogFactory.getLog(FastaWriter.class);

	public static final int DEFAULT_LINE_LENGTH = 60;

	BufferedWriter writer;
	String filename = "<unknown>";
	int lineLength;
	
	// If this instance was instantiated with a writer rather than
	// a file, this is false and this object is not responsible for
	// closing the writer
	boolean ownsWriter = true;
	
	/** Opens a fasta file for writing.  Uses the default line length.
	 * @param file Filename of the file to open
	 * @param append true if sequences should be appended to this file, false if it should be replaced. 
	 * */
	public FastaWriter(String file, boolean append) {
		this(new File(file), append, DEFAULT_LINE_LENGTH);
	}

	/** Opens a fasta file for writing.
	 * @param file Filename of the file to open
	 * @param append true if sequences should be appended to this file, false if it should be replaced. 
	 * @param lineLength Allows the length of a sequence line to be set. 
	 * */
	public FastaWriter(String file, boolean append, int lineLength) {
		this(new File(file), append, lineLength);
	}

	/** Opens a fasta file for writing.  Uses the default line length
	 * @param file Filename of the file to open
	 * @param append true if sequences should be appended to this file, false if it should be replaced. 
	 * */
	public FastaWriter(File file, boolean append) {
		this(file, append, DEFAULT_LINE_LENGTH);
	}

	/** Opens a fasta file for writing.
	 * @param file Filename of the file to open
	 * @param append true if sequences should be appended to this file, false if it should be replaced. 
	 * @param lineLength Allows the length of a sequence line to be set. 
	 * */
	public FastaWriter(File file, boolean append, int lineLength) {
		try {
			setupWriter(new FileWriter(file, append), lineLength);
			// fasta files don't require a terminal newline (are terminal newlines allowed?)
			// but initial newlines are not allowed
			if (append && file.length() > 0) {
				writer.newLine();
			}
		} catch (FileNotFoundException ex) {
			throw new ConfigException("Not able to write fasta file: " + file, ex);
		} catch (IOException ex) {
			throw new ConfigException("Error writing fasta file: " + file, ex);
		}
		filename = file.getAbsolutePath();
	}

	public FastaWriter(Writer w, int lineLength) {
		setupWriter(w, lineLength);
		ownsWriter = false;
	}

	protected void setupWriter(Writer w, int lineLength) {
		this.lineLength = lineLength;
		writer = new BufferedWriter(w);
	}
	
	/** Writes a sequence to the file. 
	 * @param header The fasta header for the sequence.  Should not include '>'.
	 * @param content The actual sequence.  It will be broken up into lines based on linelength as it is written out.
	 * */
	public void writeSeq(String header, String content) {
		try {
			writer.write(">");
			writer.write(header);
			writer.newLine();
			int i = 0;
			log.debug("Writing sequence of length: " + content.length() + " to " + filename);
			while (i < content.length()) {
				CharSequence line = content.subSequence(i, Math.min(lineLength+i, content.length()) );
				writer.write(line.toString());
				writer.newLine();
				i = i + lineLength;
			}
			writer.flush();
		}
		catch(IOException ex) {
			throw new ConfigException("Not able to write fasta file: " + filename, ex);
		}
	}

	/** Writes a portion of the sequence to the file. 
	 * @param header The fasta header for the sequence.  Should not include '>'.
	 * @param content The actual sequence.  It will be broken up into lines based on linelength as it is written out.
	 * @param start 1-based index of the first base to write.  If null, writing will start from the beginning of the sequence. 
	 * @param stop 1-based index of the last base to write.  If null, writing will continue until the end of the sequence.
	 * */
	public void writeSeq(String header, String string, Integer start, Integer stop) {
		CharSequence sequence = string;
		// trim the sequence to the specified coordinates
		// note that we don't subtract 1 from the stop coordinate,
		// because of the way String.substring works
		if (start != null && stop != null)
			sequence = sequence.subSequence(start.intValue() - 1, stop.intValue());
		else if (start != null)
			sequence = sequence.subSequence(start.intValue() - 1, sequence.length());
		else if (stop != null)
			sequence = sequence.subSequence(0, stop.intValue());
		writeSeq(header, sequence.toString());
	}

	/** Writes a sequence to the file, which may have been previously
	 * read from a different file.
	 */
	public void writeSeq(FastaSequence seq) {
		writeSeq(seq.getHeader(), seq.getSequence());
	}
	
	/** Closes the file.  Once it is closed, no more sequences may be written.
	 * Open a new FastaWriter to append more sequences into the file.
	 */
	public void close() {
		FileUtil.safeClose(writer);
	}
	
	public String getFilename() {
		return filename;
	}

	@Override
	protected void finalize() {
		if (ownsWriter) {
			close();
		}
	}
}
