/*
 * $Id:$
 *
 * $Log$
 */
/***********************************************************************
*              						               *
* MEME                                                                 *
* Copyright 2007, The University of Queensland 			       *
* Author: Timothy L. Bailey                                            *
*                                                                      *
***********************************************************************/

#include <libxslt/xslt.h>
#include <libxslt/xsltInternals.h>
#include <libxslt/transform.h>
#include <libxslt/xsltutils.h>


/**********************************************************************
print_meme_file_html

Print MEME results in HTML format.
Format XML as HTML using a stylesheet and an XSLT
**********************************************************************/
extern void print_meme_file_html(
    char* stylesheet_file_path,   /* path to MEME XSL stylesheet IN */
    char* input_file_path,        /* path to XML input file IN */
    char* output_file_path        /* path to HTML output file IN */
) {
  xsltStylesheetPtr stylesheet = NULL;
  xmlDocPtr input_doc = NULL;
  xmlDocPtr output_doc = NULL;
  const int PERFORM_ENTITY_SUBST = 1;
  xmlSubstituteEntitiesDefault(PERFORM_ENTITY_SUBST);
  xmlLoadExtDtdDefaultValue = 0;
  stylesheet = xsltParseStylesheetFile((const xmlChar *) stylesheet_file_path);
  input_doc = xmlParseFile(input_file_path);
  output_doc = xsltApplyStylesheet(stylesheet, input_doc, NULL);
  xsltSaveResultToFilename(output_file_path, output_doc, stylesheet, 0);
  xmlFreeDoc(output_doc);
  xmlFreeDoc(input_doc);
  xsltFreeStylesheet(stylesheet);

} /* print_meme_file_html */

