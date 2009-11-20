package calhoun.util;

/** Exception type for consistency check errors.  It is unchecked to simplify error handling in other Java code.
 * It is a wrapper for problems that arise from failed consistency checks.It is a subclass of ConfigException.
 * You should assume error messages will be viewed by the user.
 */
public class CheckException extends RuntimeException {

	private static final long serialVersionUID = 4225964190101139621L;

	public CheckException() {
	}

	public CheckException(String arg0) {
		super(arg0);
	}

	public CheckException(String arg0, Throwable arg1) {
		super(arg0);
		if(arg1.getCause() != null)
			initCause(arg1.getCause());
		else
			initCause(arg1);
	}

	public CheckException(Throwable arg) {
		if(arg.getCause() != null)
			initCause(arg.getCause());
		else
			initCause(arg);
	}
    
    /**
     * Gets the plain message unembellished with stack traces and causes.
     */
    public String getPlainMessage() {
        String message = super.getMessage();
        // Strip off the bean toString when it is added to the message
        int index = message.indexOf(": calhoun.");
        if ( index != -1 )
            message = message.substring(0, index);
        return message;
    }
    
	@Override
	public String getMessage() {
		String ret = super.getMessage();
		Throwable cause = getCause();
		if(cause != null) {
			ret += " caused by: "+cause.getMessage()+" - "+cause.toString();
		}
		return ret;
	}
}
