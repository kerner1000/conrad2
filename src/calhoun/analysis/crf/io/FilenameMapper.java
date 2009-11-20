package calhoun.analysis.crf.io;

import java.io.File;

/** An interface used to map filenames.  Given a source filename, converts it into a destination filename. */
public interface FilenameMapper {

	/** maps a filename from the source to a destination.
	 * 
	 * @param source the input filename to use as the base for the mapping
	 * @return the destination filename
	 */
	File mapFilename(File source);
}
