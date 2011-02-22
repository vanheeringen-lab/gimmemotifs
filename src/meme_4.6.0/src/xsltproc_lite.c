/**********************************************************************
  xsltproc_lite

  This program transforms an XML file by applying an XSLT styleshett
  and write the output to stdout

**********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <libxslt/xslt.h>
#include <libxslt/xsltInternals.h>
#include <libxslt/transform.h>
#include <libxslt/xsltutils.h>

int main(int argc, char* argv[]) {

  char* usage = "xsltpro_lite <xslt filename> <xml input filename> <output filename>\n";

  if (argc != 4) {
    fprintf(stderr, "%s", usage);
    exit(0);
  }

  char* stylesheet_file_path = argv[1]; // Path to XSLT stylesheet
  char* input_file_path = argv[2];      // Path to XML input file IN
  char* output_file_path = argv[3];      // Path to XML input file IN

  xsltStylesheetPtr stylesheet = NULL;
  xmlDocPtr input_doc = NULL;
  xmlDocPtr output_doc = NULL;
  const int PERFORM_ENTITY_SUBST = 1;
  xmlSubstituteEntitiesDefault(PERFORM_ENTITY_SUBST);
  xmlLoadExtDtdDefaultValue = 0;
  stylesheet = xsltParseStylesheetFile((xmlChar *) stylesheet_file_path);
  if (! stylesheet) {
    fprintf(stderr, "Couldn't parse stylesheet file %s\n", stylesheet_file_path);
    exit(1);
  }
  input_doc = xmlParseFile(input_file_path);
  if (! stylesheet) {
    fprintf(stderr, "Couldn't parse input file %s\n", input_file_path);
    exit(1);
  }
  output_doc = xsltApplyStylesheet(stylesheet, input_doc, NULL);
  if (! output_doc) {
    fprintf(stderr, "Couldn't apply stylesheet %s to input file %s\n",
      stylesheet_file_path, input_file_path);
    exit(1);
  }
  xsltSaveResultToFilename(output_file_path, output_doc, stylesheet, 0);
  xmlFreeDoc(output_doc);
  xmlFreeDoc(input_doc);
  xsltFreeStylesheet(stylesheet);

  exit(0);

}; 
