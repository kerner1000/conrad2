package calhoun.seq;

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Iterator;
import java.util.zip.GZIPInputStream;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import calhoun.util.DataException;

public class QualityIterator implements Iterator {
	private static final Log log = LogFactory.getLog(QualityIterator.class);

	File file;
	BufferedReader reader;
	
	String lastReadHeader;
	byte [] lastReadQualitySequence;
	
	String lastReadLine;
	byte [] lastReadQuality;
	boolean endOfFileReached;

	public QualityIterator(String fileName, boolean isGzipped) throws IOException {
		this(new File(fileName), isGzipped);
	}
	
	public QualityIterator(String fileName) throws IOException {
		this(new File(fileName), false);
	}
	
	public QualityIterator(BufferedReader reader) {
		commonInit(reader);
	}
	
	protected void commonInit(BufferedReader reader) {
		endOfFileReached = false;
		this.reader = reader;
		readNext();
	}
	
	public QualityIterator(File file, boolean isGzipped) throws IOException {
		// save the file for error reporting
		this.file = file;
		
		InputStream is = new FileInputStream(file);
		
		if(isGzipped) {
			is = new GZIPInputStream(is);
		}
		
		// create a BufferedReader to read the file
		reader = new BufferedReader(new InputStreamReader(is));
		commonInit(reader);
	}

	/** Returns the next set of qualities in this file. */
	public Object next() {
		return nextQuality();
	}

	/** Returns the next set of qualities in this file.  Provides a detailed return type. */
	public FastaSequence nextQuality() {
		SimpleFastaSequence qs = new SimpleFastaSequence();
		addNextQuality(qs);
		return qs;
	}
	
	/** Adds the next sequence of qualities in the file to the given FastaSequence.  If the FastaSequence already has a sequence set, verifies that the quality matches the sequence (header and length). */
	public void addNextQuality(SimpleFastaSequence input) {
		if(endOfFileReached)
			throw new DataException("Attempted to read past end of file "+file.getAbsolutePath());

		// If the header is set on the input sequence, verify that it is the same.
		if(input.header == null)
			input.setHeader(lastReadHeader);
		else if(!input.header.equals(lastReadHeader)) {
			throw new DataException("Quality file header "+lastReadHeader+" doesn't match sequence header "+input.header);
		}

		// If the sequence is set on the input sequence, verify that it is the same length.
		if(input.hasSequence() && (input.getLength() != lastReadQualitySequence.length)) {
			throw new DataException("Quality file length "+lastReadQualitySequence.length+" doesn't match sequence length "+input.getLength()+" for "+input.header);
		}
			
		input.setQuality(lastReadQualitySequence);
		
		readNext();
	}
	
	public boolean hasNext() {
		return !endOfFileReached;
	}
	
	protected void readNext() {
		String line;
		try {
		do {
			if(lastReadLine == null)
				line = reader.readLine();
			else
				line = lastReadLine;
			
			if(line == null) {
				endOfFileReached = true;
				return;
			}
			
			line = line.trim();
			if(!line.startsWith(">")) {
				throw new DataException("Expected \">\" when reading next quality sequence");
			}
		} while(line.length() == 0);

		lastReadHeader = line.substring(1).trim();
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		
		while(true) {
			line = reader.readLine();
			if(line == null) {
				lastReadLine = null;
				break;
			}
			
			line = line.trim();
			if(line.startsWith(">")) {
				lastReadLine = line;
				break;
			}
			
			if(line.length() > 0) {
				String numbers [] = line.split(" ");
				for(int i=0;i<numbers.length;i++) {
					int val = Integer.parseInt(numbers[i]);
					os.write(val);
				}
			}
		}
		
		lastReadQualitySequence = os.toByteArray();
		} catch(IOException e) {
			throw new DataException("Exception reading "+file.getPath(), e);
		}
	}
	
	public void remove() {
		throw new UnsupportedOperationException();
	}
}
