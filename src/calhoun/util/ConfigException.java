package calhoun.util;

/** Exception type for system configuration exceptions.  It is unchecked to simplify error handling in other Java code.
 * It is a wrapper for problems that generally will not be handled in code.  Files not found, setup missing, etc.
 * You should assume error messages will be viewed by the user.
 * 
 * This should only be thrown for errors which should be fixable by the user without examining or changing code. 
 */
public class ConfigException extends RuntimeException{

	private static final long serialVersionUID = -5371694423739150838L;

	public ConfigException(String arg0) {
		super(arg0);
	}

	public ConfigException(String arg0, Throwable arg1) {
		//super(arg0, arg1);
		super(arg0);
		if(arg1.getCause() != null)
			initCause(arg1.getCause());
		else
			initCause(arg1);
	}

	public ConfigException(Throwable arg) {
		//super(arg0);
		if(arg.getCause() != null)
			initCause(arg.getCause());
		else
			initCause(arg);
	}
}
