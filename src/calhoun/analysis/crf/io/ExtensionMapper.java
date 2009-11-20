package calhoun.analysis.crf.io;

import java.io.File;

import calhoun.util.Assert;
import calhoun.util.FileUtil;

/** an implementation of {@link FilenameMapper} that changes file extensions or add a new file extension onto the source file.
 * If the <code>append<code> property is set, the new extension will be added to the end of the existing filename.  Otherwise,
 * the existing file extension is replaced with the new extension.
 */
public class ExtensionMapper implements FilenameMapper {

	boolean append;
	String extension;

	/** returns the value of the append property.  If true, the <code>extension<code> will be appended to the source filename.  
	 * Otherwise, the existing extension will replaced with the new one.  Default is false.
	 * @return the value of the append property.
	 */
	public boolean isAppend() {
		return append;
	}

	/** set the value of the append property.
	 * @param append set to true if the extension should be appended to create the new filename.
	 */
	public void setAppend(boolean append) {
		this.append = append;
	}
	
	/** the new extension to use in the mapped file.  Must be set before {link #mapFilename} is called.
	 * @return the extension to use in the mapped file.
	 */
	public String getExtension() {
		return extension;
	}
	
	/** sets the extension.
	 * @param extension the extension to use in the mapped file.
	 */
	public void setExtension(String extension) {
		this.extension = extension;
	}

	/** Maps the source file name to a new name by altering the file extension.  The new <code>extension</code>
	 * will either be added to the end of the existing filename or will replace the existing extension based on the
	 * value of the <code>append</code> property.
	 * @see calhoun.analysis.crf.io.FilenameMapper#mapFilename(java.io.File)
	 */
	public File mapFilename(File source) {
		Assert.a(extension != null, "Extension must be set befure mapFilename can be called.");
		String dest;
		if(append) {
			dest = source.getPath()+extension;
		}
		else {
			String[] baseAndExtension = FileUtil.getBaseAndExtension(source);
			String directory = "";
			if(source.getParent() != null) {
				directory = source.getParent()+File.separator;
			}
			dest = directory+baseAndExtension[0]+extension;
		}
		return new File(dest);
	}
}
