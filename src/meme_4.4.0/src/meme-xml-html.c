/*
 * $Id:$
 * 
 * $Log$
 */

/************************************************************************
*	Copyright							*
*	(2007) The University of Queensland				*
*	All Rights Reserved.						*
*	Author: Timothy L. Bailey					*
*									*
*	Permission to use, copy, modify, and distribute any part of 	*
*	this software for educational, research and non-profit purposes,*
*	without fee, and without a written agreement is hereby granted, *
*	provided that the above copyright notice, this paragraph and 	*
*	the following three paragraphs appear in all copies.		*
*									*
*	Those desiring to incorporate this software into commercial 	*
*	products or use for commercial purposes should contact the 	*
*	Technology Transfer Office, University of California, San Diego,*
*	9500 Gilman Drive, La Jolla, California, 92093-0910, 		*
*	Ph: (619) 534 5815.						*
*									*
*	IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO 	*
*	ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR 	*
*	CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF 	*
*	THE USE OF THIS SOFTWARE, EVEN IF THE UNIVERSITY OF CALIFORNIA 	*
*	HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 		*
*									*
*	THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE *
*	UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE 		*
*	MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.  *
*	THE UNIVERSITY OF CALIFORNIA MAKES NO REPRESENTATIONS AND 	*
*	EXTENDS NO WARRANTIES OF ANY KIND, EITHER EXPRESSED OR IMPLIED, *
*	INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 	*
*	MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT 	*
*	THE USE OF THE MATERIAL WILL NOT INFRINGE ANY PATENT, 		*
*	TRADEMARK OR OTHER RIGHTS.  					*
************************************************************************/

#define DEFINE_GLOBALS
#include "meme-print-html.h"
#include "macros.h"
#include "general.h"

// Convert a MEME XML file into an HTML file.
int main(int argc, char **argv) {

  char *xml_file_name = NULL;
  char *xsl_file_name = NULL;
  char *html_file_name = NULL;

  DO_STANDARD_COMMAND_LINE(1,
    USAGE([options]);
    DATA_OPTN(1, xml, <xml>, \tname of xml file, xml_file_name = _OPTION_);
    DATA_OPTN(1, xsl, <xsl>, \tname of xsl file, xsl_file_name = _OPTION_);
    DATA_OPTN(1, html, <html>, \tname of html file, html_file_name = _OPTION_);
    USAGE(\n\tConvert MEME XML to HTML using the given style sheet.);
    USAGE(\n\tCopyright);
    USAGE(\t(2007) The University of Queensland);
    USAGE(\tAll Rights Reserved.);
    USAGE(\tAuthor: Timothy L. Bailey);
  );

  if (!xml_file_name) { 
    fprintf(stderr, "You must specify an XML file.\n");
    exit(1);
  }
  if (!xsl_file_name) { 
    fprintf(stderr, "You must specify an XSL file.\n");
    exit(1);
  }
  if (!html_file_name) { 
    fprintf(stderr, "You must specify an HTML file.\n");
    exit(1);
  }
  
  print_meme_file_html(xsl_file_name, xml_file_name, html_file_name);

  return(0);
} // main
