package calhoun.util;


/** Exception type for general errors.  It is unchecked to simplify error handling in the rest of our java code.
 * It should be used for unrecoverable errors that are probably not due to a system configuration.  Users should
 * be expected to find a developer when one of these occurs.  In general these will be bugs in the code or corner 
 * cases you have not implemented.  If there is a whole operation or method that is completely not working, just 
 * throw UnsupportedOperationException.  
 */
public final class ErrorException extends RuntimeException{
	
	private static final long serialVersionUID = -7560596902065799827L;

	public ErrorException(String arg0) {
		super(arg0);
	}

	public ErrorException(String arg0, Throwable arg1) {
		//super(arg0, arg1);
		super(arg0);
		if(arg1.getCause() != null)
			initCause(arg1.getCause());
		else
			initCause(arg1);
	}

	public ErrorException(Throwable arg) {
		//super(arg0);
		if(arg.getCause() != null)
			initCause(arg.getCause());
		else
			initCause(arg);
	}

	@Override
	public String getMessage() {
		StringBuffer ret = new StringBuffer(); 
		ret.append(super.getMessage());
		Throwable cause = getCause();
		if(cause != null) {
			ret.append(" caused by: ").append(cause.getMessage()).append(" - ").append(cause.toString());
			/*if(cause.getMessage() == null) {
				Writer w = new StringWriter();
				cause.printStackTrace(new PrintWriter(w));
				ret.append(w);
			}*/
		}
		return ret.toString();
	}
}
