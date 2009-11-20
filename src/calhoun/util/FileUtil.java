package calhoun.util;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.io.Reader;
import java.io.Writer;
import java.nio.channels.FileChannel;
import java.nio.channels.FileLock;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

/** Utility functions for file management
 */
public class FileUtil {
	private static final Log log = LogFactory.getLog(FileUtil.class);

	/** This class should not be instantiated. */
	private FileUtil() {
	}

	public static final String readFile(File filename) throws IOException {
		return readReader(new FileReader(filename));
	}

	public static final void appendSeparator(File in) throws Exception {
		FileOutputStream f = new FileOutputStream(in, true);
		f.write(java.io.File.separator.getBytes());
	}
	public static final void copyFile(File in, File out, boolean append) throws Exception {
		FileChannel sourceChannel = new FileInputStream(in).getChannel();
		FileChannel destinationChannel = new FileOutputStream(out, append).getChannel();
		sourceChannel.transferTo(0, sourceChannel.size(), destinationChannel);
		sourceChannel.close();
		destinationChannel.close();
	}
	public static final void copyFile(File in, File out) throws Exception {
		copyFile(in, out, false);
	}
	
	public static final void renameFile(File oldFile, File newFile) throws Exception {
		newFile.getParentFile().mkdirs();
		if(!oldFile.renameTo(newFile)) {
			copyFile(oldFile, newFile);
			oldFile.delete();
		}
	}
	
	public static final BufferedWriter safeOpen(final String file) {
		try {
			return file == null ? null :
				new BufferedWriter(new FileWriter(file));
		}
		catch(IOException ex) {
			throw new RuntimeException(ex);
		}
	}
	
	public static final void safeWrite(final Writer w, final String s) {
		try {
			w.write(s);
		}
		catch(IOException ex) {
			throw new RuntimeException(ex);
		}
	}
	
	public static final void unzipFile(File file, File toDirectory) throws Exception {
	    toDirectory.mkdirs();
	    ZipInputStream zis = new ZipInputStream(new FileInputStream(file));
	    ZipEntry entry;
	    while ( ( entry = zis.getNextEntry() ) != null ) {
	        if ( entry.isDirectory() )
	            new File(toDirectory, entry.getName()).mkdirs();
	        else {
	            FileOutputStream fos = new FileOutputStream(new File(toDirectory, entry.getName()));
	            byte[] bytes = new byte[8096];
	            int read;
	            while ( ( read = zis.read(bytes) ) != -1 ) {
	                fos.write(bytes, 0, read);
	            }
	        }
	    }
	}
	
	/**
	 * Checks to see if the stream is in GZIP format, assuming that
	 * it is open to the start of a file or other resource.  Note that
	 * this advances the stream, caller should mark and reset around this method.
	 */
	public static final boolean isGzipStream(InputStream s) {
		try {
			byte[] b = {0,0};
			int magic = 0;
			if (s.read(b) == 2) {
				magic = ((b[1]&0xff) << 8 | (b[0]&0xff));
			}
			return magic == GZIPInputStream.GZIP_MAGIC;
		}
		catch (IOException e) {
			throw new ErrorException(e);
		}		
	}
	
	public static final boolean isGzipFile(File file) {
		try {
			InputStream s = new FileInputStream(file);
			boolean result = isGzipStream(s);
			s.close();
			return result;
		}
		catch (IOException e) {
			throw new ErrorException(e);
		}
	}
	
	public static final boolean isGzipFile(String filename) {
		return isGzipFile(new File(filename));
	}
	
	public static final String[] getBaseAndExtension(File file) {
		String name = file.getName();
		int period = name.lastIndexOf('.');
		if(period == -1) {
			return new String[] { name, null };
		}
		else {
			return new String[] { name.substring(0, period), name.substring(period+1) };
		}
	}
	
	public static File makeTempCopy(File in) throws Exception{
		int index = 0;
		String pathPrefix = in.getAbsolutePath();
		File f;
		do {
			f = new File(pathPrefix + "temp_"+index);
			index ++;
		} while(f.exists());
		
		copyFile(in, f);
		
		return f;
	}
	
	public static final String readInputStream(InputStream is) throws IOException {
		BufferedReader reader = new BufferedReader(new InputStreamReader(is));
		return readReader(reader);
	}
	
	public static final String readReader(Reader reader) throws IOException {
		char buf[] = new char[4096];
		StringBuffer ret = new StringBuffer();
		int size = 0;
		do {
			size = reader.read(buf);
			if (size != -1)
				ret.append(buf, 0, size);
		} while (size != -1);
		return ret.toString();
	}

	public static final byte[] readFileAsBytes(String filename) throws IOException {
		File f = new File(filename);
		byte[] data = new byte[(int) f.length()];
		InputStream is = new FileInputStream(f);
		is.read(data);
		is.close();
		return data;
	}

	/** Reads a files into a String.
	 * 
	 * @param filename  Name of the file to read in
	 * @return String containing the contents of the file
	 */
	public static final String readFile(String filename) throws IOException {
		return readFile(new File(filename));
	}

	public static String[][] readFlatFile(String fileName) throws IOException {
		BufferedReader in = new BufferedReader(new FileReader(fileName));
		ArrayList<String[]> ret = new ArrayList<String[]>();
		String line;
		while ((line = in.readLine()) != null) {
//			System.out.println(line);
			if(line.trim().charAt(0) == '#') {
				continue;
			}
			ret.add(line.split("\t"));
		}
		return (String[][]) ret.toArray(new String[ret.size()][]);
	}
	
	public static String[][] readFlatFileWithComments(String fileName) throws IOException {
		BufferedReader in = new BufferedReader(new FileReader(fileName));
		ArrayList<String[]> ret = new ArrayList<String[]>();
		String line;
		while ((line = in.readLine()) != null) {
			ret.add(line.split("\t"));
		}
		return (String[][]) ret.toArray(new String[ret.size()][]);
	}	
	
	public static double[] readDoublesFromSingleTabbedLine(String fileName) throws IOException {
		String[][] s = FileUtil.readFlatFile(fileName);
		Assert.a(s.length == 1);
		int len = s[0].length;
		double[] ret = new double[len];
		for (int j=0; j<len; j++) { ret[j] = Double.parseDouble(s[0][j]) ; }
		return ret;
	}

	public static String[] readLines(String fileName) throws IOException {
		BufferedReader in = new BufferedReader(new FileReader(fileName));
		ArrayList ret = new ArrayList();
		String line;
		while ((line = in.readLine()) != null) {
			ret.add(line);
		}
		return (String[]) ret.toArray(new String[ret.size()]);
	}

	/** Writes a string out to a file.
	 * 
	 * @param filename  Name of the file to read in
	 * @param data The string to write to the file
	 */
	public static final void writeFile(File filename, String data) throws IOException {
		FileWriter fw = new FileWriter(filename);
		PrintWriter logWriter = new PrintWriter(fw);
		logWriter.print(data);
		fw.close();
	}

	/** Serializes a Java object out to a file.
	 * 
	 * @param filename  Name of the file to write
	 * @param data The object to serialize.  Must be Serializable
	 */
	public static final void writeObject(String filename, Object data) throws IOException {
		OutputStream fw = new BufferedOutputStream(new FileOutputStream(filename));
		ObjectOutputStream objectWriter = new ObjectOutputStream(fw);
		objectWriter.writeObject(data);
		objectWriter.close();
	}

	/** Read a serialized Java object from a file.
	 * 
	 * @param filename  Name of the file to read in
	 * @return The object which was deserialized.
	 */
	public static final Object readObject(String filename) throws IOException, ClassNotFoundException {
		FileInputStream fr = new FileInputStream(filename);
		ObjectInputStream objectIn = new ObjectInputStream(fr);
		Object object = objectIn.readObject();
		objectIn.close();
		return object;
	}

	/** Compares records in two flat files.  Just checks that every line in one file occurs in the other.  Reports differences.
	 * Like fileDiff, except that it is insensitive to differences in ordering.
	 */
	public static final String fileRecordCompare(String expected, String actual, String maskingRegEx) throws IOException {
		BufferedReader expectedReader = new BufferedReader(new FileReader(expected));
		
		Map records = new HashMap();
		while (true) {
			String s = expectedReader.readLine();
			if (s == null)
				break;
			if(maskingRegEx != null)
				s = s.replaceAll(maskingRegEx, "[!MASKED!]");
			Number count = (Number)records.get(s);
			if (count == null) records.put(s,new Integer(1));
			else records.put(s,new Integer(count.intValue()+1));
		}
		expectedReader.close();

		StringBuffer ret = new StringBuffer();
		BufferedReader actualReader = new BufferedReader(new FileReader(actual));
		while (true) {
			String s = actualReader.readLine();
			if (s == null)
				break;
			if(maskingRegEx != null)
				s = s.replaceAll(maskingRegEx, "[!MASKED!]");
			if (records.containsKey(s)) {
				Number count = (Number)records.get(s);
				if (count.intValue() == 1) records.remove(s);
				else records.put(s,new Integer(count.intValue()-1));
			} else {
				ret.append("+ " + s + "\n");
			}
		}
		for (Iterator iter = records.keySet().iterator(); iter.hasNext();) {
			String element = (String) iter.next();
			ret.append("- " + element + "\n");
		}
		if (ret.length() == 0) {
			return null;
		} else {
			return ret.toString();
		}
	}

	/** Convience form of fileRecordCompare which does no masking
	 */
	public static final String fileRecordCompare(String expected, String actual) throws IOException {
		return fileRecordCompare(expected, actual, null);
	}
	
	/** Writes a string out to a file.
	 * 
	 * @param filename  Name of the file to read in
	 * @param data The string ot write to the file
	 */
	public static final void writeFile(String filename, String data) throws IOException {
		writeFile(new File(filename), data);
	}

	/** Adds a file or directory name to an existing directory name.  Handles duplicate or missing slashes
	 * 
	 * @param parent Name of the parent directory
	 * @param child  Name of the child directory
	 * @return The concatenated path.
	 */
	public static final String appendPath(String parent, String child) {
		// Use forward separator and not the system default.  Forward works on windows.
		String separator = "/";
		StringBuffer buf = new StringBuffer(parent);
		if (!parent.endsWith(separator))
			buf.append(separator);
		if (child.startsWith(separator))
			buf.append(child.substring(1));
		else
			buf.append(child);
		return buf.toString();
	}

	public static final void safeClose(Writer f) {
		try {
			if (f != null)
				f.close();
		} catch (IOException ex) {
			throw new ErrorException("Error closing writer: " + f, ex);
		}
	}

	public static final void safeClose(Reader f) {
		try {
			if (f != null)
				f.close();
		} catch (IOException ex) {
			throw new ErrorException("Error closing reader: " + f, ex);
		}
	}

	public static final void safeClose(OutputStream f) {
		try {
			if (f != null)
				f.close();
		} catch (IOException ex) {
			throw new ErrorException("Error closing outputStream: " + f, ex);
		}
	}

	public static final void safeClose(InputStream f) {
		try {
			if (f != null)
				f.close();
		} catch (IOException ex) {
			throw new ErrorException("Error closing inputStream: " + f, ex);
		}
	}

	public static final void safeClose(RandomAccessFile f) {
		try {
			if (f != null)
				f.close();
		} catch (IOException ex) {
			throw new ErrorException("Error closing randomAccessFile: " + f, ex);
		}
	}

	/**
	 * Locks a file.
	 *
	 * This lock will work across JVMs but not necessarily across threads
	 * within a single JVM. You must have write permission to lock a file.
	 */
	public static final FileLock lockFile(RandomAccessFile raf) {
		// given a RandomAccessFile object, gets the unique channel associated
		// with it and locks it
		FileChannel ch = raf.getChannel();

		// after the lock has been acquired, continue to use the raf object to write
		// data to the File since RAF's implement different data output methods than the
		// FileChannel class. The lock object will behave essentially as a mutex: if 
		// you have it, you've got the lock. If you don't, wait for it. 

		// return the lock object - technically not necessary as the lock will 
		// be released when the raf is closed or when FileUtil.safeRelease gets called
		// but it might be useful to have around
		return doLock(ch);
	}

	/**
	 * Locks a file.
	 *
	 * This lock will work across JVMs but not necessarily across threads
	 * within a single JVM. You must have write permission to lock a file.
	 */
	public static final FileLock lockFile(FileOutputStream fos) {
		FileChannel ch = fos.getChannel();
		return doLock(ch);
	}
	
	private static final FileLock doLock(FileChannel ch) {
		FileLock lock;
		try {
			lock = ch.lock();
		} catch (IOException ex) {
			log.warn("File locking not supported on this platform. Continuing.");
			log.warn("If the file " + ch + " is in use by another process, " +
			  "silent corruption may occur.");
			lock = null;
		} catch (Exception ex) {
			throw new ErrorException("Error trying to acquire lock on file. ", ex);
		}
		return lock;
	}
	
	public static final void safeRelease(FileLock lock) {
		try {
			if (lock != null)
				lock.release();
		} catch (IOException ex) {
			throw new ErrorException("Error releasing lock: " + lock, ex);
		}
	}
	
	/** Determines the size of a newline character in the given text file.  Will return 1 or 2 depending on the
	 * platform the file was created on.  
	 * @param file the file object to examine.
	 * @return the size of a newline in this file
	 * @throws IOException
	 */
	public static byte determineNewlineSize(File file) throws IOException {
		InputStream stream = new FileInputStream(file);
		try {
			byte[] data = new byte[1000];

			int count;
			do {
				count = stream.read(data);
				log.debug("byte read in " + file.getName() +" : " + count );
				for (int i = 0; i < count; ++i) {
					if (data[i] == 10) {
						log.debug("Found Newline");
						return 1;
					} else if (data[i] == 13) {
						log.debug("Found CarriageReturn");
						// Handle the edge case where 13 is the last character in the buffer 
						int nextChar = i + 1 < count ? data[i + 1] : stream.read();
						return (byte) ((nextChar == 10) ? 2 : 1);
					}
				}
			} while (count != -1);
			// if File does not contain a newLine, we can safely assume that the size of a newline is 1.
			return 1;
			//throw new ConfigException("File did not contain a newline, could not determine type.");
		} finally {
			stream.close();
		}
	}
}
