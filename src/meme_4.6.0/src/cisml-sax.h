/********************************************************************
 * FILE: cisml-sax.h
 * AUTHOR: James Johnson
 * CREATE DATE: 7 September 2010
 * PROJECT: MEME suite
 * COPYRIGHT: 2010, UQ
 ********************************************************************/
#ifndef CISML_SAX_H
#define CISML_SAX_H

#include <libxml/parser.h>

typedef struct CISML_CALLBACKS CISML_CALLBACKS_T;

struct CISML_CALLBACKS {
  void (*start_cisml)(void*);                               // <status>
  void (*start_cis_element_search)(void*);                  // <status>
  void (*handle_program_name)(void*, char*);                // <status> <program name>
  void (*start_parameters)(void*);                          // <status>
  void (*handle_pattern_file)(void*, char*);                // <status> <pattern file>
  void (*handle_sequence_file)(void*, char*);               // <status> <sequence file>
  void (*handle_background_seq_file)(void*, char*);         // <status> <background sequence file>
  void (*handle_pattern_pvalue_cutoff)(void*, double);      // <status> <pattern pvalue cutoff>
  void (*handle_sequence_pvalue_cutoff)(void*, double);     // <status> <sequence pvalue cutoff>
  void (*handle_site_pvalue_cutoff)(void*, double);         // <status> <site pvalue cutoff>
  void (*handle_sequence_filtering)(void*, int, char*);     // <status> <filtering used?> [<filtering program>]    
  void (*end_parameters)(void*);                            // <status>
  void (*start_multi_pattern_scan)(void*, double*, double*);// <status> [<pvalue>] [<score>]
  // read the start of a pattern
  // <status> <accession> <name> [<db>] [<lsId>] [<pvalue>] [<score>]
  void (*start_pattern)(void*, char*, char*, char*, char*, double*, double*);
  // read the start of a scanned sequence
  // <status> <accession> <name> [<db>] [<lsId>] [<score>] [<pvalue>] [<length>]
  void (*start_scanned_sequence)(void*, char*, char*, char*, char*, double*, double*, long*);
  // read the matched element
  // <status> <start> <stop> [<score>] [<pvalue>] [<clusterId>]
  void (*start_matched_element)(void*, long, long, double*, double*, char*);
  // read the matched sequence for the immediately preceeding call to handle_matched_element
  // <status> <sequence>
  void (*handle_sequence)(void*, char*);
  void (*end_matched_element)(void*);
  void (*end_scanned_sequence)(void*);                      // <status>
  void (*end_pattern)(void*);                               // <status>
  void (*end_multi_pattern_scan)(void*);                    // <status>
  void (*end_cis_element_search)(void*);                    // <status>
  void (*end_cisml)(void*);                                 // <status>
  //for extension to CISML pass on SAX parser information
  void (*start_unknown)(void*, const xmlChar*, const xmlChar**);
  void (*end_unknown)(void*, const xmlChar*);
  void (*characters_unknown)(void*, const xmlChar*, int len);
};

int parse_cisml(CISML_CALLBACKS_T *callbacks, void *state, const char *filename);

#endif

