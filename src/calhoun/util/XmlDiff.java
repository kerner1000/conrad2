/* The Broad Institute SOFTWARE COPYRIGHT NOTICE AGREEMENT
This software and its documentation are copyright 2003 by the Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
This software is supplied without any warranty or guaranteed support whatsoever. Neither the Broad Institute nor MIT can be responsible for its use, misuse, or functionality. 
*/
package calhoun.util;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.regex.Pattern;

import org.apache.commons.lang.StringUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.dom4j.Attribute;
import org.dom4j.Document;
import org.dom4j.Element;

/** Utility class for diffing XML
 */
public class XmlDiff {
	private static final Log log = LogFactory.getLog(XmlDiff.class);

	/** Should not be constructed.  All static methods.
	 * 
	 */
	private XmlDiff() {
		super();
	}

	/** Compares two XML files.  the first is the template, which can contain '*' as a wildcard attribute or element value.  Returns a String describing the differences or null if the files match. 
	 * The comparison is an XML aware comparison that ignores whitespace differences and attribute order and also allows wildcard values.  
		Just put '*' as the value for any element content or attribute value and any plugin output will match.  It also ignores ordering of 
		duplicated child elements.  Therefore if <BlastRun> has 50 <BlastAlignment> child objects, the order does not have to match in the 
		2 documents.*/
	public static String compareFiles(String templateName, String docName) {
		log.debug("Comparing files: "+templateName+" and "+docName);
		Document template = XmlUtil.parseFile(templateName);
		Document doc = XmlUtil.parseFile(docName);
		List result = new ArrayList();

		if(template.getRootElement().getName() != doc.getRootElement().getName())
			result.add("Root nodes don't match.  Expected: '"+template.getRootElement().getName()+"', Received: "+doc.getRootElement().getName());
		else
			result = compareXmlElements(template.getRootElement(), doc.getRootElement());

		if(result.size() == 0)
			return null;
		else
			return "Files differ: "+docName+" doesn't match "+templateName+"\n"+StringUtils.join(result.iterator(), '\n');
	}

	static List compareXmlElements(Element template, Element doc) {
		List result = new ArrayList();
		
		log.debug("Comparing elements: "+template.getUniquePath()+" and "+doc.getUniquePath());

		Set attributes = new HashSet(doc.attributes());

		// Compare attributes
		Iterator it = template.attributes().iterator();
		while(it.hasNext()) {
			Attribute attribute = (Attribute) it.next(); 
			Attribute docAttribute = doc.attribute(attribute.getQName());
			if(docAttribute == null)
				result.add(template.getUniquePath()+": Expected: "+attribute.getQualifiedName()+"='"+attribute.getValue()+"', Received: No attribute");
			else {
				String comp = compareValues(attribute.getValue(), docAttribute.getValue());
				if(comp != null)
					result.add(attribute.getUniquePath()+": "+comp);
				attributes.remove(docAttribute);
			}
		}
		it = attributes.iterator();
		while(it.hasNext()) {
			Attribute attribute = (Attribute) it.next(); 
			result.add(template.getUniquePath()+": Expected: No Attribute, Received: "+attribute.getQualifiedName()+"='"+attribute.getValue()+"'");
		}

		// Compare content
		if(template.isTextOnly()) {
			if(!doc.isTextOnly())
				result.add(template.getUniquePath()+": Expected: '"+template.getText()+"', Received: Non-text content");
			else {
				String comp = compareValues(template.getText(), doc.getText());
				if(comp != null)
					result.add(template.getUniquePath()+": "+comp);
			}
		}


		Set matchedElements = new HashSet();
		// Compare children
		it = template.elements().iterator();
		while(it.hasNext()) {
			Element child = (Element) it.next(); 
			List docElements = new ArrayList(doc.elements(child.getQName()));
			docElements.removeAll(matchedElements);
			if(docElements.size() == 0) {
				result.add(template.getUniquePath()+": Expected: <"+child.getQualifiedName()+">, Received: No child element");
			}
			else {
				// Loop through children, saving the best matching element.  Stop if we get a perfect match.
				int fewestProbs=999999999;
				Iterator docIit = docElements.iterator();
				Element bestElement = null;
				List bestResults = null;
				while(docIit.hasNext()) {
					Element docElement = (Element) docIit.next(); 
					List childResult = compareXmlElements(child, docElement);
					if(childResult.size() < fewestProbs) {
						fewestProbs = childResult.size();
						bestElement = docElement;
						bestResults = childResult;
					}
					if(fewestProbs == 0)
						break;
				}
				// Now take the best element as the match
				matchedElements.add(bestElement);
				result.addAll(bestResults);
				if(docElements.size() > 1)
					log.debug("Matched: "+child.getUniquePath()+" with "+bestElement.getUniquePath());
			}
		}				
		Set unmatchedElements = new HashSet(doc.elements());
		unmatchedElements.removeAll(matchedElements);
		it = unmatchedElements.iterator();
		while(it.hasNext()) {
			Element element = (Element) it.next(); 
			result.add(template.getUniquePath()+": Expected: No Child Element, Received: <"+element.getQualifiedName()+">");
		}
	
		return result;
	}

	static String compareValues(String templateValue, String docValue) {
		if (templateValue.equals("*")) {
			return null;
		}
		else if (templateValue.startsWith("REGEX:")) {
			String regex = templateValue.substring("REGEX:".length());
			if (!Pattern.compile(regex).matcher(docValue).matches()) {
				return "Expected: '"+templateValue+"', Received: '"+docValue+"'";				
			}
		}
		else if(!templateValue.equals(docValue)) {
			return "Expected: '"+templateValue+"', Received: '"+docValue+"'";
		}
		return null;
	}
}
