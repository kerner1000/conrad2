package calhoun.util;

/** Exception type for data errors.  It is unchecked to simplify error handling in other Java code.
 * It is a wrapper for problems that are generally attributable to bad data in the database.  It is a subclass of ConfigException.
 * You should assume error messages will be viewed by the user.
 */
public class DataException extends ConfigException {

	private static final long serialVersionUID = 1993656669514976055L;

	public DataException(String arg0) {
		super(arg0);
	}

	public DataException(String arg0, Throwable arg1) {
		super(arg0, arg1);
	}

	public DataException(Throwable arg0) {
		super(arg0);
	}
}
