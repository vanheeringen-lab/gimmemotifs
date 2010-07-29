#ifndef BEADSTRING_XML_H
#define BEADSTRING_XML_H
#include "mhmm-state.h"
#include "beadstring.h"
#include "cisml.h"
#include <stdio.h>
void print_beadstring_xml_header(FILE *out);
void print_beadstring_xml_alphabet(FILE *out, MHMM_T *hmm);
void print_beadstring_xml_model(FILE *out, OPTIONS_T *options, MHMM_T  *hmm);
void print_beadstring_xml_cisml(FILE *out, MHMM_T *hmm, CISML_T *cisml);
void print_beadstring_xml_trailer(FILE *out);
#endif
