package calhoun.util;


import java.io.PrintWriter;
import java.io.StringWriter;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.logging.Formatter;
import java.util.logging.LogRecord;

/** Simpler formatter than the default Java SimpleFormatter.  Only displays time and file.
 */
public class SimplerFormatter extends Formatter {

	static SimpleDateFormat f = new SimpleDateFormat("HH:mm:ss");
	static String blanks="       ";

	/** Just prints a time followed by message followed by class function line
	 */
	@Override
	public String format(LogRecord rec) {
		String l= rec.getLevel().getName();
		String err = "";
		Throwable t = rec.getThrown();
		if(t != null) {
			err = t.getMessage() + ": "+printStackTrace(t);
		}
		return f.format(new Date(rec.getMillis()))+blanks.substring(0, 7-l.length())
		+l+" "+rec.getLoggerName()+"."+rec.getSourceMethodName()+" - "+rec.getMessage()+"\n"+err;
	}

	private String printStackTrace(Throwable t) {
		StringWriter sw = new StringWriter();
		t.printStackTrace(new PrintWriter(sw, true));
		return sw.toString().replaceAll("\\p{Blank}*at org.python.core.*\n", "");
	}
}
