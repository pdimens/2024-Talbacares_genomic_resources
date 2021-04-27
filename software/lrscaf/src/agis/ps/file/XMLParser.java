/*
*File: agis.ps.util.XMLParser.java
*User: mqin
*Email: mqin@ymail.com
*Date: 2016年1月22日
*/
package agis.ps.file;

import java.io.File;
import java.io.IOException;
//import java.util.regex.Matcher;
//import java.util.regex.Pattern;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import agis.ps.util.Parameter;

public class XMLParser {
	private static Logger logger = LoggerFactory.getLogger(XMLParser.class);

	public Parameter parseXML(String file) {
		Parameter para = null;
		try {
			DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
			DocumentBuilder builder = factory.newDocumentBuilder();
			Document doc = builder.parse(new File(file));
			para = new Parameter();
			// root element
			Element rootElm = doc.getDocumentElement();
//			logger.debug(this.getClass().getName() + "\t" + rootElm.getTagName());

			NodeList nodes = rootElm.getChildNodes();
			if (nodes == null || nodes.getLength() == 0) {
				logger.error(this.getClass().getName() + "\t" + "The XML file do not contain any elements!");
				return null;
			}

			// traverse child input elements
			nodes = rootElm.getElementsByTagName("input");
			if (nodes == null || nodes.getLength() == 0) {
				logger.error(this.getClass().getName() + "\t" + "The XML file do not contain input elements!");
				return null;
			} else {
				Node input = nodes.item(0);
				nodes = input.getChildNodes();
				for (int i = 0; i < nodes.getLength(); i++) {
					Node node = nodes.item(i);
					if (node.getNodeType() == Node.ELEMENT_NODE) {
						String nodeName = ((Element) node).getNodeName();
						if (nodeName.equalsIgnoreCase("contig")) {
							para.setCntFile(node.getTextContent().trim());
						} else if (nodeName.equalsIgnoreCase("m5")) {
							para.setAlgFile(node.getTextContent().trim());
							para.setType("m5");
						} else if (nodeName.equalsIgnoreCase("m4")){
							para.setAlgFile(node.getTextContent().trim());
							para.setType("m4");
						}else if (nodeName.equalsIgnoreCase("sam")) {
							para.setAlgFile(node.getTextContent().trim());
							para.setType("sam");
						} else if(nodeName.equalsIgnoreCase("mm")){
							para.setAlgFile(node.getTextContent().trim());
							para.setType("mm");
						} else if(nodeName.equalsIgnoreCase("pbread")){
							para.setPbFile(node.getTextContent().trim());
						} else {
							logger.info(this.getClass().getName() + "\t" + "The XML file contain illeagle elements!");
							return null;
						}
					}
				}
			}

			// traverse child output elements
			nodes = rootElm.getElementsByTagName("output");
			if (nodes == null || nodes.getLength() == 0) {
				logger.info(this.getClass().getName() + "\t" + "The XML file do not contain output elements, it will use "
						+ System.getProperty("user.dir") + "!");
				para.setOutFolder(System.getProperty("user.dir"));
			} else {
				Node node = nodes.item(0);
				para.setOutFolder(node.getTextContent().trim());
			}

			// traverse child paras elements
			nodes = rootElm.getElementsByTagName("paras");
			if (nodes == null || nodes.getLength() == 0) {
				logger.info(this.getClass().getName() + "\t" + "The XML file do not contain output elements, it will use default value!");
				para.setMinContLen(200);
				para.setMinPBLen(5000);
				para.setMinOLLen(160);
				para.setMinOLRatio(0.8);
				para.setMaxOHLen(300);
				para.setMaxOHRatio(0.1);
				para.setMaxEndLen(300);
				para.setMaxEndRatio(0.1);
				para.setMinSupLinks(1);
				para.setMaxSupLinks(50);
				if(para.getType() == "mm")
				{
					para.setIdentity(0.1);
				} else
				{
					para.setIdentity(0.8);
				}
				para.setUseOLLink(true);
				para.setRatio(0.2);
				para.setRepMask(true);
				para.setGapFilling(false);
				para.setTipLength(1500);
				para.setIqrTime(1.5);
				para.setThreads(4);
			} else {
				Node parasNode = nodes.item(0);
				nodes = parasNode.getChildNodes();
				for (int i = 0; i < nodes.getLength(); i++) {
					Node node = nodes.item(i);
					if (node.getNodeType() != Node.ELEMENT_NODE)
						continue;
					String nodeName = node.getNodeName();
					if (nodeName.equalsIgnoreCase("min_contig_length")) {
						int value = Integer.valueOf(node.getTextContent().trim());
						if(value < 200)
						{
							para.setMinContLen(200);
							logger.warn("The minimum contig's length should be large than 200 bp!");
							logger.warn("Mandatorily set to 200 bp!");
						} else
						{
							para.setMinContLen(value);
						}
					} else if (nodeName.equalsIgnoreCase("min_pacbio_length")) {
						para.setMinPBLen(Integer.valueOf(node.getTextContent().trim()));
					} else if (nodeName.equalsIgnoreCase("min_overlap_length")) {
						para.setMinOLLen(Integer.valueOf(node.getTextContent().trim()));
					} else if (nodeName.equalsIgnoreCase("min_overlap_ratio")) {
						para.setMinOLRatio(Double.valueOf(node.getTextContent().trim()));
					} else if (nodeName.equalsIgnoreCase("min_supported_links")) {
						para.setMinSupLinks(Integer.valueOf(node.getTextContent().trim()));
					} else if (nodeName.equalsIgnoreCase("max_supported_links")) {
						para.setMaxSupLinks(Integer.valueOf(node.getTextContent().trim()));
					} else if (nodeName.equalsIgnoreCase("max_overhang_length")) {
						para.setMaxOHLen(Integer.valueOf(node.getTextContent().trim()));
					} else if (nodeName.equalsIgnoreCase("max_overhang_ratio")) {
						para.setMaxOHRatio(Double.valueOf(node.getTextContent().trim()));
					} else if (nodeName.equalsIgnoreCase("max_end_length")) {
						para.setMaxEndLen(Integer.valueOf(node.getTextContent().trim()));
					} else if (nodeName.equalsIgnoreCase("max_end_ratio")) {
						para.setMaxEndRatio(Double.valueOf(node.getTextContent().trim()));
					} else if (nodeName.equalsIgnoreCase("identity")) {
						// para.setIdentity(Double.valueOf(node.getTextContent().trim()));
						double ident = Double.valueOf(node.getTextContent().trim());
						if(para.getType().equalsIgnoreCase("mm") && ident >= 0.3)
						{
							logger.warn("The identity is setted to " + ident + ".");
							logger.warn("The identity for Minimap mapper should be not large than 0.3!");
						} 
						para.setIdentity(ident);
					} else if (nodeName.equalsIgnoreCase("use_overlap_link")) {
						String temp = node.getTextContent().trim();
						if(temp.startsWith("t") || temp.startsWith("T"))
							para.setUseOLLink(true);
						else 
							para.setUseOLLink(false);
					} else if (nodeName.equalsIgnoreCase("ratio")){
						para.setRatio(Double.valueOf(node.getTextContent().trim()));
					} else if (nodeName.equalsIgnoreCase("repeat_mask")){
						String temp = node.getTextContent().trim();
						if(temp.startsWith("t") || temp.startsWith("T"))
							para.setRepMask(true);
						else 
							para.setRepMask(false);
					} else if(nodeName.equalsIgnoreCase("gap_filling")){
						String temp = node.getTextContent().trim();
						if(temp.startsWith("t") || temp.startsWith("T"))
							para.setGapFilling(true);
						else 
							para.setGapFilling(false);
					} else if(nodeName.equalsIgnoreCase("tip_length")){
						para.setTipLength(Integer.valueOf(node.getTextContent().trim()));
					} else if(nodeName.equalsIgnoreCase("iqr_time")){
						para.setIqrTime(Double.valueOf(node.getTextContent().trim()));
					} else if(nodeName.equalsIgnoreCase("mmcm")){
						para.setMmcm(Integer.valueOf(node.getTextContent().trim()));;
					} else if(nodeName.equals("process")){
						para.setThreads(Integer.valueOf(node.getTextContent().trim()));
					} else {
						logger.info(this.getClass().getName() + "\t" + "The para element contain illeage item " + nodeName + ". it will be omitted!");
					}
				}
			}

		} catch (ParserConfigurationException e) {
			logger.debug("Error: ", e);
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
		} catch (SAXException e) {
			logger.debug("Error: ", e);
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
		} catch (IOException e) {
			logger.debug("Error: ", e);
			logger.error(this.getClass().getName() + "\t" + e.getMessage());
		}
		return para;
	}
}
