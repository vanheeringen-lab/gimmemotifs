/**************************************************************************
 * FILE: xml-io.h
 * CREATE DATE: 8/3/2007
 * AUTHOR: William Stafford Noble, Tim Bailey, and Charles Grant
 * PROJECT: MEME
 * COPYRIGHT: 2007, University of Washington
 * DESCRIPTION: Utility functions for reading XML files.
 **************************************************************************/
#ifndef XML_UTIL_H
#define XML_UTIL_H
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>
#include <libxslt/xslt.h>
#include <libxslt/xsltInternals.h>
#include <libxslt/transform.h>
#include <libxslt/xsltutils.h>
#include "array.h"

/***********************************************************************
 * Look up a node set using an XPath string.
 * Caller is responsible for freeing the xmlXPathObj.
 ***********************************************************************/
xmlXPathObjectPtr xpath_query(
  xmlXPathContextPtr xpath_ctxt, 
  char* path
);

/***********************************************************************
 * Check if property from an XML node exists.
 * Caller is responsible for freeing the xmlChar string.
 ***********************************************************************/
BOOLEAN_T check_xml_node_property(
  xmlNodePtr node, 
  char* property_name
);

/***********************************************************************
 * Read a property from an XML node.
 * Caller is responsible for freeing the xmlChar string.
 ***********************************************************************/
xmlChar* read_xml_node_property(
  xmlNodePtr node, 
  char* property_name
);

/***********************************************************************
 * Read and set the alphabet from a parsed XML document.
 ***********************************************************************/
void read_alphabet_from_xml(xmlXPathContextPtr xpath_ctx);

/***********************************************************************
 * Read the background letter frequencies from XML.
 * Caller is responsible for freeing the returned array.
 ***********************************************************************/
ARRAY_T *read_bg_freqs_from_xml (xmlXPathContextPtr xpath_ctx);

/**********************************************************************
print_xml_to_file_using_stylesheet

Print the contents of an XML Doc to a file, applying an 
XSLT stylesheet.
**********************************************************************/
BOOLEAN_T print_xml_filename_to_file_using_stylesheet(
    char* input_file_path,      /* path to XML input file IN */
    char* stylesheet_file_path, /* path to MEME XSL stylesheet IN */
    FILE* output_file           /* path to HTML output file IN */
);

/**********************************************************************
print_xml_file_to_filename_using_stylesheet

Print the contents of an XML file to another file applying an 
XSLT stylesheet.

Returns TRUE if successful, FALSE otherwise.
**********************************************************************/
BOOLEAN_T print_xml_filename_to_filename_using_stylesheet(
    char* input_file_path,       /* path to XML input file IN */
    char* stylesheet_file_path,  /* path to MEME XSL stylesheet IN */
    char* output_file_path       /* path to HTML output file IN */
);

#endif
