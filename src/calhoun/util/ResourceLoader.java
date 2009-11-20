package calhoun.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.StringWriter;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

/**
 * The ResourceLoader loads for resources by trying each of the following mechanisms:
 * <ul>
 * <li>Treating the resourceID as a file name and opening it in the file system
 * <li>Loading a class resource with the specified resourceID
 * <li>Loading a system resource with the specified resourceID
 * </ul>  
 */
public class ResourceLoader {

	private static final Log log = LogFactory.getLog(ResourceLoader.class);

	private static ResourceFinder[] resourceFinders = { new FileSystem(), new ClassResource(), new SystemResource() }; 

	/**
	 * Opens an InputStream by trying each of the defined ResourceFinders.
	 */	
	public static InputStream openInputStream(Class cls, String resourceID)
	{	
		log.debug("Attempting to load resource "+ resourceID);
		resourceID = resourceID.replace('\\', '/');
		for ( int i = 0; i < resourceFinders.length; ++i ) {
			log.debug("\tTrying resource finder "+ resourceFinders[i]);
			InputStream is = resourceFinders[i].openResource(cls, resourceID);
			if ( is != null ) {
				log.debug("\tResource found");
				return is;			
			}
		}
		throw new ConfigException("Resource '" + resourceID + "' not found for " + cls);
	}

	public static String openTextResource(Class cls, String resourceId) {
		try {
			BufferedReader reader = new BufferedReader(new InputStreamReader(openInputStream(cls, resourceId)));
			StringWriter writer = new StringWriter();
			PrintWriter print = new PrintWriter(writer);
			String line;
			while ( ( line = reader.readLine() ) != null ) {
				print.println(line);
			}
			print.flush();
			return writer.toString();
		}
		catch (IOException x) {
			throw new ErrorException(x);
		}
	}

	public interface ResourceFinder
	{
		/**
		 * @return null if the resource cannot be found
		 */
		public InputStream openResource(Class cls, String resourceID);
	}
	
	static class FileSystem
		implements ResourceFinder
	{
		public InputStream openResource(Class cls, String resourceID)
		{
			File file = new File(resourceID);
			if ( !file.exists() || !file.isFile() )
				return null;
			try
			{
				return new FileInputStream(file);
			}
			catch (FileNotFoundException x)
			{
				log.warn("Unable to open existing file as resource", x);
				return null;
			}
		}		
	}
	
	static class ClassResource
		implements ResourceFinder
	{
		public InputStream openResource(Class cls, String resourceID)
		{
			return cls.getResourceAsStream(resourceID);
		}		
	}

	static class SystemResource
		implements ResourceFinder
	{
		public InputStream openResource(Class cls, String resourceID)
		{
			return cls.getClassLoader().getResourceAsStream(resourceID);
		}		
	}
}
