/*
 * The Broad Institute
 * 
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * 
 * This software and its documentation are copyright 2004 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 * 
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse, or functionality. 
 * 
 * Created on Jul 2, 2004
 */
package calhoun.seq;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Iterator;
import java.util.NoSuchElementException;

import calhoun.util.ConfigException;
import calhoun.util.ErrorException;

/**
 * This class implements iterator and allows a fasta file to be 
 * iterated through without having to read anything.  Can also take a second file containing quality.
 * 
 * You can now use new FastaReader().iterator() instead.
 */
public class FastaIterator implements Iterator<FastaSequence> {
	
	public static final String HEADER_START = ">";
	
	private File file;	
	private BufferedReader reader;	
	private String	nextHeader;
	private QualityIterator	qualityIt;

	public FastaIterator(String file) throws IOException {
		this(new File(file), null);
	}
	
	public FastaIterator(File file) throws IOException {
		this(file, null);
	}
	
	public FastaIterator(String sequence, String quality) throws IOException {
		this(new File(sequence), quality == null ? null : new File(quality));
	}

	/** Creates an iterator over a sequence file and also possibly a quality file */
	public FastaIterator(File sequence, File quality) throws IOException {
		
		// save the file for error reporting
		this.file = sequence;
		if(!file.exists())
			throw new IOException("The fasta file \""+file+"\" does not exist");
		
		if(file.length() == 0) {
			nextHeader = null;
			return;
		}
	
		if(quality != null) {
			qualityIt = new QualityIterator(quality, false);
		}
		
		// create a BufferedReader to read the file
		reader = new BufferedReader(new FileReader(file));
		
		// read in the first line
		String nextHeaderLine;
		do {
			nextHeaderLine = reader.readLine();
		} while(nextHeaderLine.trim().length() == 0);
		
		// validate the line
		if (nextHeaderLine.startsWith(HEADER_START) == false)
			throw new ConfigException(file + " is not a valid FASTA file");
			
		// clean up the header line
		nextHeader = nextHeaderLine.substring(1).trim();
	}
	
	/** Throws an UnsupportedOperationException. */
	public void remove() {
		throw new UnsupportedOperationException("can't remove items from a FastaIterator");
	}

	public boolean hasNext() {		
		return nextHeader != null;
	}
	
	public FastaSequence next() {
		
		if (nextHeader == null)
			throw new NoSuchElementException();
		
		ReaderSequence seq = new ReaderSequence(file, nextHeader);
		
		try {
			nextHeader = seq.loadSequence(reader);
		}
		catch (Exception e) {
			throw new RuntimeException("error reading " + this.file, e);
		}
		
		try {
			if(nextHeader == null)
				reader.close();
		} catch (IOException ex) {
			throw new ErrorException("Error closing "+file.getAbsolutePath(), ex);
		}
		
		if(qualityIt != null) {
			qualityIt.addNextQuality(seq);
		}
		return seq;
	}
	
	/** A minimal sequence class that slurps its sequence in immediately. */

	public static class ReaderSequence extends SimpleFastaSequence {
		private static final long serialVersionUID = -637754971987954577L;
		File file;
		
		public ReaderSequence (File file, String header) {
			this.file = file;
			setHeader(header);
		}
		
		public File getFile() {
			return file;
		}
		
		public String loadSequence(BufferedReader reader) throws IOException {
			
			// read the sequence
				
				String line;
				StringBuffer buffer = new StringBuffer();
				
				while (true) {
					
					// get the next line from the reader
							line = reader.readLine();
						
					// see if the line starts a new sequence
						if (line == null || line.startsWith(HEADER_START))
							break;
						else
							buffer.append(line.trim());
				}
			
			setSequence(buffer.toString());
			
			// extract the header from the header line
				if (line != null)
					line = line.substring(1).trim();
			
			return line;
		}
	}
}
