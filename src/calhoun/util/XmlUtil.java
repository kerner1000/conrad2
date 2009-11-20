package calhoun.util;

import java.io.BufferedInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.StringWriter;
import java.util.zip.GZIPInputStream;

import javax.xml.transform.Templates;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.stream.StreamResult;
import javax.xml.transform.stream.StreamSource;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.dom4j.Document;
import org.dom4j.DocumentException;
import org.dom4j.DocumentHelper;
import org.dom4j.Element;
import org.dom4j.Node;
import org.dom4j.io.OutputFormat;
import org.dom4j.io.SAXReader;
import org.dom4j.io.XMLWriter;

/** Utility functions to simplify XML.  Wrapper around DOM4j.
 */
public final class XmlUtil {
	@SuppressWarnings("unused")
	private static final Log log = LogFactory.getLog(XmlUtil.class);

	/** This class should not be instantiated. */
	private XmlUtil() {
	}

	public static void transform(String input, String output, String stylesheet) {
		try {
			TransformerFactory tFactory = TransformerFactory.newInstance();
			StreamSource source = new StreamSource(ResourceLoader.openInputStream(XmlUtil.class, stylesheet));
			Templates templates = tFactory.newTemplates(source);
			Transformer transformer = templates.newTransformer();
			StreamSource sourceXML = new StreamSource(input);
			StreamResult result = new StreamResult(output);
			
			transformer.transform(sourceXML, result);	
		}
		catch(Exception e) {
			throw new ConfigException(e);
		}
	}
	
	public static Document newDocument() {
		return DocumentHelper.createDocument();
	}
		
	public static Element newElement(String name) {
		return DocumentHelper.createElement(name);
	}
	
	/** Creates a dom4j Document from a file.  The filename given can be on the classpath or loaded from the filesystem.
	 * Converts exceptions into runtime ConfigExceptions.  
	 */
	public static Document parseFile(String filename) {
		try {
			SAXReader xmlReader = new SAXReader();
			InputStream s = new BufferedInputStream(ResourceLoader.openInputStream(XmlUtil.class,filename));
			s.mark(100);
			boolean isGzip = FileUtil.isGzipStream(s);
			s.reset();
			if (isGzip) s = new GZIPInputStream(s);
			return xmlReader.read(s);
		} catch (IOException e) {
			throw new ErrorException(e);
		} catch (DocumentException ex) {
			throw new ConfigException("Invalid XML file: ", ex);
		}
	}

	/** Creates a dom4j Document from a text string.
	 * Converts exceptions into runtime ErrorExceptions.  
	 */
	public static Document parseString(String data) {
		try {
			return DocumentHelper.parseText(data);
		} catch (DocumentException ex) {
			throw new ErrorException("Invalid XML data: ", ex);
		}
	}

	/** Pretty prints a document and returns the result as a String */
	public static String prettyPrint(Document doc) {
		return prettyPrint((Node) doc);
	}

	/** Pretty prints a dom4j Node and returns the result as a String */
	public static String prettyPrint(Node node) {
		StringWriter sw = new StringWriter();
		try {
			OutputFormat outformat = OutputFormat.createPrettyPrint();
			XMLWriter writer = new XMLWriter(sw, outformat);
			writer.write(node);
			writer.flush();
		} catch (Exception ex) {
			throw new ErrorException("Error writing XML output ", ex);
		}
		return sw.toString();
	}

	/** Pretty prints a document to a given OutputStream */
	public static void prettyPrint(Document doc, OutputStream out) {
		try {
			OutputFormat outformat = OutputFormat.createPrettyPrint();
			XMLWriter writer = new XMLWriter(out, outformat);
			writer.write(doc);
			writer.flush();
		} catch (Exception ex) {
			throw new ErrorException("Error writing XML output ", ex);
		}
	}

	/** Pretty prints a document to a given filename */
	public static void prettyPrint(Document doc, String filename) {
		try {
			OutputStream out = new FileOutputStream(filename);
			prettyPrint(doc, out);
			out.close();
		} catch (IOException ex) {
			throw new ConfigException("Error saving XML to file ", ex);
		}
	}

}
