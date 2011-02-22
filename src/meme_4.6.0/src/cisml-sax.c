/********************************************************************
 * FILE: cisml-sax.c
 * AUTHOR: James Johnson
 * CREATE DATE: 7 September 2010
 * PROJECT: MEME suite
 * COPYRIGHT: 2010, UQ
 ********************************************************************/

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdio.h>

#include "binary-search.h"
#include "cisml-sax.h"
#include "utils.h"

/*debugging macros {{{*/
extern VERBOSE_T verbosity;

#define DEBUG_MSG(debug_level, debug_msg) { \
  if (verbosity >= debug_level) { \
    fprintf(stderr, debug_msg); \
  } \
}

#define DEBUG_FMT(debug_level, debug_msg_format, ...) { \
  if (verbosity >= debug_level) { \
    fprintf(stderr, debug_msg_format, __VA_ARGS__); \
  } \
}
/*}}}*/

/*parser state {{{*/

#define PS_ERROR 0
#define PS_START 1
#define PS_IN_CIS_ELEMENT_SEARCH 2
#define PS_IN_PROGRAM_NAME 3
#define PS_IN_PARAMETERS 4
#define PS_IN_PATTERN_FILE 5
#define PS_IN_SEQUENCE_FILE 6
#define PS_IN_BACKGROUND_SEQ_FILE 7
#define PS_IN_PATTERN_PVALUE_CUTOFF 8
#define PS_IN_SEQUENCE_PVALUE_CUTOFF 9
#define PS_IN_SITE_PVALUE_CUTOFF 10
#define PS_IN_SEQUENCE_FILTERING 11
#define PS_IN_MULTI_PATTERN_SCAN 12
#define PS_IN_PATTERN 13
#define PS_IN_SCANNED_SEQUENCE 14
#define PS_IN_MATCHED_ELEMENT 15
#define PS_IN_SEQUENCE 16
#define PS_END 17

char const * const state_names[] = {
  "ERROR",
  "START",
  "IN_CIS_ELEMENT_SEARCH",
  "IN_PROGRAM_NAME",
  "IN_PARAMETERS",
  "IN_PATTERN_FILE",
  "IN_SEQUENCE_FILE",
  "IN_BACKGROUND_SEQ_FILE",
  "IN_PATTERN_PVALUE_CUTOFF",
  "IN_SEQUENCE_PVALUE_CUTOFF",
  "IN_SITE_PVALUE_CUTOFF",
  "IN_SEQUENCE_FILTERING",
  "IN_MULTI_PATTERN_SCAN",
  "IN_PATTERN",
  "IN_SCANNED_SEQUENCE",
  "IN_MATCHED_ELEMENT",
  "IN_SEQUENCE",
  "END"
};

#define MP_UNDECIDED 0
#define MP_SINGLE_PATTERN 1
#define MP_MULTI_PATTERN 2

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Keeps track of the parser state to ensure the CISML is valid
// and provides storage for callbacks
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

typedef struct CHARBUF CHARBUF_T;
typedef struct CISML_PARSER CISML_PARSER_T;

struct CHARBUF {
  int (*accept) (char);
  char *buffer;
  int pos;
  int size;
};

struct CISML_PARSER {
  CISML_CALLBACKS_T *callbacks;
  void *invoker_state;
  int state;
  int multi;
  int udepth;
  CHARBUF_T characters;
};


/*}}}*/

/* parsing helper functions {{{*/

typedef struct multi MULTI_T;
struct multi {
  int count;
  char **options;
  int *outputs;
  int *target;
};

int compare_pstrings(const void *v1, const void *v2) {
  char *str1, *str2;
  str1 = *((char**)v1);
  str2 = *((char**)v2);
  return strcmp(str1, str2);
}

int ld_str(char *value, void *data) {
  *((char**)data) = value;
  return 0;
}

int ld_char(char *value, void *data) {
  int i;
  char c;
  for (i = 0; value[i] != '\0'; ++i) {
    if (!isspace(value[i])) break;
  }
  //check that we found something
  if (value[i] == '\0') {
    return -1;
  }
  c = value[i];
  //check that it was the only thing to be found
  for (++i; value[i] != '\0'; ++i) {
    if (!isspace(value[i])) return -2;
  }
  //set the value
  *((char*)data) = c;
  return 0;
}

int ld_int(char *value, void *data) {
  long parsed_value;
  char *end_ptr;
  errno = 0;
  parsed_value = strtol(value, &end_ptr, 10);
  if (end_ptr == value) return -1;
  if (errno) return errno;
  if (parsed_value > INT_MAX || parsed_value < INT_MIN) return ERANGE;
  *((int*)data) = (int)parsed_value;
  return 0;
}

int ld_long(char *value, void *data) {
  long parsed_value;
  char *end_ptr;
  errno = 0;
  parsed_value = strtol(value, &end_ptr, 10);
  if (end_ptr == value) return -1;
  if (errno) return errno;
  *((long*)data) = parsed_value;
  return 0;
}

int ld_double(char *value, void *data) {
  double parsed_value;
  char *end_ptr;
  errno = 0;
  parsed_value = strtod(value, &end_ptr);
  if (end_ptr == value) return -1;
  if (errno) return errno;
  *((double*)data) = parsed_value;
  return 0;
}

int ld_pvalue(char *value, void *data) {
  double parsed_value;
  char *end_ptr;
  errno = 0;
  parsed_value = strtod(value, &end_ptr);
  if (end_ptr == value) return -1;
  if (errno) return errno;
  if (parsed_value < 0 || parsed_value > 1) return -1;
  *((double*)data) = parsed_value;
  return 0;
}

int ld_multi(char *value, void *data) {
  int found;
  MULTI_T *attr;
  attr = (MULTI_T*)data;
  found = binary_search(&value, attr->options, attr->count, sizeof(char*), compare_pstrings);
  if (found >= 0) {
    *(attr->target) = attr->outputs[found];
    return 0;
  } else {
    return -1;
  }
}

void parse_attributes(CISML_PARSER_T *ps, char *tag, const xmlChar **attrs, 
    int count, char **names, int (**parsers)(char*, void*), void **parser_data, 
    BOOLEAN_T *required, BOOLEAN_T *done) {
  int i, found;
  BOOLEAN_T done_data[count];
  if (done == NULL) done = done_data;
  for (i = 0; i < count; ++i) done[i] = FALSE;
  for (i = 0; attrs[i] != NULL; i += 2) {
    found = binary_search(attrs+i, names, count, sizeof(char*), compare_pstrings);
    if (found >= 0) {
      int (*parser)(char*, void*);
      void* data;
      if (done[found]) {
        ps->state = PS_ERROR;
        DEBUG_FMT(HIGH_VERBOSE, "CISML parser error: Duplicate attribute %s::%s.\n", tag, names[found]);
        continue;
      }
      done[found] = TRUE;
      parser = parsers[found];
      data = parser_data[found];
      if (parser((char*)attrs[i+1], data)) {
        ps->state = PS_ERROR;
        DEBUG_FMT(HIGH_VERBOSE, "CISML parser error: Bad value \"%s\" for attribute %s::%s.\n", (char*)attrs[i+1], tag, names[found]);
        continue;
      }
    }
  }
  for (i = 0; i < count; ++i) {
    if (required[i] && !done[i]) {
      ps->state = PS_ERROR;
      DEBUG_FMT(HIGH_VERBOSE, "CISML parser error: Missing required attribute %s::%s.\n", tag, names[i]);
    }
  }
}
/*}}}*/

void start_ele_cis_element_search(CISML_PARSER_T *ps, const xmlChar **attrs) {
  //we're not interested in any of the attributes
  if (ps->callbacks->start_cis_element_search) {
    ps->callbacks->start_cis_element_search(ps->invoker_state);
  }
}

void end_ele_program_name(CISML_PARSER_T *ps) {
  if (ps->callbacks->handle_program_name) {
    ps->callbacks->handle_program_name(ps->invoker_state, ps->characters.buffer);
  }
}

void start_ele_parameters(CISML_PARSER_T *ps, const xmlChar **attrs) {
  if (ps->callbacks->start_parameters) {
    ps->callbacks->start_parameters(ps->invoker_state);
  }
}

void end_ele_pattern_file(CISML_PARSER_T *ps) {
  if (ps->callbacks->handle_pattern_file) {
    ps->callbacks->handle_pattern_file(ps->invoker_state, ps->characters.buffer);
  }
}

void end_ele_sequence_file(CISML_PARSER_T *ps) {
  if (ps->callbacks->handle_sequence_file) {
    ps->callbacks->handle_sequence_file(ps->invoker_state, ps->characters.buffer);
  }
}

void end_ele_background_seq_file(CISML_PARSER_T *ps) {
  if (ps->callbacks->handle_background_seq_file) {
    ps->callbacks->handle_background_seq_file(ps->invoker_state, ps->characters.buffer);
  }
}

void end_ele_pattern_pvalue_cutoff(CISML_PARSER_T *ps) {
  double pvalue;
  char *end_ptr;
  pvalue = strtod(ps->characters.buffer, &end_ptr);
  if (ps->characters.buffer == end_ptr || *end_ptr != '\0' || pvalue < 0 || pvalue > 1) {
    ps->state = PS_ERROR;
    DEBUG_MSG(HIGH_VERBOSE, "CISML parser error: pattern-pvalue-cutoff must be number in range 0 to 1 inclusive");
    return;
  }
  if (ps->callbacks->handle_pattern_pvalue_cutoff) {
    ps->callbacks->handle_pattern_pvalue_cutoff(ps->invoker_state, pvalue);
  }
}

void end_ele_sequence_pvalue_cutoff(CISML_PARSER_T *ps) {
  double pvalue;
  char *end_ptr;
  pvalue = strtod(ps->characters.buffer, &end_ptr);
  if (ps->characters.buffer == end_ptr || *end_ptr != '\0' || pvalue < 0 || pvalue > 1) {
    ps->state = PS_ERROR;
    DEBUG_MSG(HIGH_VERBOSE, "CISML parser error: sequence-pvalue-cutoff must be number in range 0 to 1 inclusive");
    return;
  }
  if (ps->callbacks->handle_sequence_pvalue_cutoff) {
    ps->callbacks->handle_sequence_pvalue_cutoff(ps->invoker_state, pvalue);
  }
}

void end_ele_site_pvalue_cutoff(CISML_PARSER_T *ps) {
  double pvalue;
  char *end_ptr;
  pvalue = strtod(ps->characters.buffer, &end_ptr);
  if (ps->characters.buffer == end_ptr || *end_ptr != '\0' || pvalue < 0 || pvalue > 1) {
    ps->state = PS_ERROR;
    DEBUG_MSG(HIGH_VERBOSE, "CISML parser error: site-pvalue-cutoff must be number in range 0 to 1 inclusive");
    return;
  }
  if (ps->callbacks->handle_site_pvalue_cutoff) {
    ps->callbacks->handle_site_pvalue_cutoff(ps->invoker_state, pvalue);
  }
}

void start_ele_sequence_filtering(CISML_PARSER_T *ps, const xmlChar **attrs) {
  int is_filtered;
  char *program = NULL;

  char* filter_options[2] = {"off", "on"};
  int filter_values[2] = {FALSE, TRUE};
  MULTI_T filter_multi = {.count = 2, .options = filter_options, .outputs = filter_values, .target = &(is_filtered)};

  char* names[2] = {"on-off", "type"};
  int (*parsers[2])(char*, void*) = {ld_multi, ld_str};
  void *data[2] = {&filter_multi, &program};
  BOOLEAN_T required[2] = {TRUE, FALSE};

  parse_attributes(ps, "sequence-filtering", attrs, 2, names, parsers, data, required, NULL);

  if (ps->callbacks->handle_sequence_filtering) {
    ps->callbacks->handle_sequence_filtering(ps->invoker_state, is_filtered, program);
  }
}

void end_ele_parameters(CISML_PARSER_T *ps) {
  if (ps->callbacks->end_parameters) {
    ps->callbacks->end_parameters(ps->invoker_state);
  }
}

void start_ele_multi_pattern_scan(CISML_PARSER_T *ps, const xmlChar **attrs) {
  if (ps->multi == MP_UNDECIDED) ps->multi = MP_MULTI_PATTERN;
  else if (ps->multi == MP_SINGLE_PATTERN) {
    DEBUG_MSG(HIGH_VERBOSE, "CISML parser error: when multi-pattern-scan is used, pattern must not be a child of cis-element-search\n");
    ps->state = PS_ERROR;
    return;
  }
  double pvalue;
  double score;

  char* names[2] = {"pvalue", "score"};
  int (*parsers[2])(char*, void*) = {ld_pvalue, ld_double};
  void *data[2] = {&pvalue, &score};
  BOOLEAN_T required[2] = {FALSE, FALSE};
  BOOLEAN_T done[2] = {FALSE, FALSE};

  parse_attributes(ps, "multi-pattern-scan", attrs, 2, names, parsers, data, required, done);

  if (ps->callbacks->start_multi_pattern_scan) {
    ps->callbacks->start_multi_pattern_scan(ps->invoker_state, (done[0] ? &pvalue : NULL), (done[1] ? &score : NULL));
  }
}

void start_ele_pattern(CISML_PARSER_T *ps, const xmlChar **attrs) {
  if (ps->multi == MP_UNDECIDED) ps->multi = MP_SINGLE_PATTERN;
  else if (ps->multi == MP_MULTI_PATTERN) {
    DEBUG_MSG(HIGH_VERBOSE, "CISML parser error: when multi-pattern-scan is used, pattern must not be a child of cis-element-search\n");
    ps->state = PS_ERROR;
    return;
  }
  char *accession = NULL;
  char *name = NULL;
  char *db = NULL;
  char *lsId = NULL;
  double score = 0;
  double pvalue = 0;

  char* names[6] = {"accession", "db", "lsId", "name", "pvalue", "score"};
  int (*parsers[6])(char*, void*) = {ld_str, ld_str, ld_str, ld_str, ld_pvalue, ld_double};
  void *data[6] = {&accession, &db, &lsId, &name, &pvalue, &score};
  BOOLEAN_T required[6] = {TRUE, FALSE, FALSE, TRUE, FALSE, FALSE};
  BOOLEAN_T done[6] = {FALSE, FALSE, FALSE, FALSE, FALSE, FALSE};

  parse_attributes(ps, "pattern", attrs, 6, names, parsers, data, required, done);

  if (ps->callbacks->start_pattern) {
    ps->callbacks->start_pattern(ps->invoker_state, accession, name, db, lsId, 
        (done[4] ? &pvalue : NULL), (done[5] ? &score : NULL));
  }
}

void start_ele_scanned_sequence(CISML_PARSER_T *ps, const xmlChar **attrs) {
  char *accession = NULL;
  char *name = NULL;
  char *db = NULL;
  char *lsId = NULL;
  double score = 0;
  double pvalue = 0;
  long length = 0;

  char* names[7] = {"accession", "db", "length", "lsId", "name", "pvalue", "score"};
  int (*parsers[7])(char*, void*) = {ld_str, ld_str, ld_long, ld_str, ld_str, ld_pvalue, ld_double};
  void *data[7] = {&accession, &db, &length, &lsId, &name, &pvalue, &score};
  BOOLEAN_T required[7] = {TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE};
  BOOLEAN_T done[7] = {FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE};

  parse_attributes(ps, "scanned-sequence", attrs, 7, names, parsers, data, required, done);

  if (ps->callbacks->start_scanned_sequence) {
    ps->callbacks->start_scanned_sequence(ps->invoker_state, accession, name, db, lsId, 
        (done[6] ? &score : NULL), (done[5] ? &pvalue : NULL), (done[2] ? &length : NULL));
  }
}

void start_ele_matched_element(CISML_PARSER_T *ps, const xmlChar **attrs) {
  long start = 0;
  long stop = 0;
  double score = 0;
  double pvalue = 0;
  char *clusterId = NULL;

  char *names[5] = {"clusterId", "pvalue", "score", "start", "stop"};
  int (*parsers[5])(char*, void*) = {ld_str, ld_pvalue, ld_double, ld_long, ld_long};
  void *data[5] = {&clusterId, &pvalue, &score, &start, &stop};
  BOOLEAN_T required[5] = {FALSE, FALSE, FALSE, TRUE, TRUE};
  BOOLEAN_T done[5] = {FALSE, FALSE, FALSE, FALSE, FALSE};

  parse_attributes(ps, "matched-element", attrs, 5, names, parsers, data, required, done);

  if (ps->callbacks->start_matched_element) {
    ps->callbacks->start_matched_element(ps->invoker_state, start, stop, (done[2] ? &score : NULL), (done[1] ? &pvalue : NULL), clusterId);
  }
}

void end_ele_sequence(CISML_PARSER_T *ps) {
  if (ps->callbacks->handle_sequence) {
    ps->callbacks->handle_sequence(ps->invoker_state, ps->characters.buffer);
  }
}

void end_ele_matched_element(CISML_PARSER_T *ps) {
  if (ps->callbacks->end_matched_element) {
    ps->callbacks->end_matched_element(ps->invoker_state);
  }
}

void end_ele_scanned_sequence(CISML_PARSER_T *ps) {
  if (ps->callbacks->end_scanned_sequence) {
    ps->callbacks->end_scanned_sequence(ps->invoker_state);
  }
}

void end_ele_pattern(CISML_PARSER_T *ps) {
  if (ps->callbacks->end_pattern) {
    ps->callbacks->end_pattern(ps->invoker_state);
  }
}

void end_ele_multi_pattern_scan(CISML_PARSER_T *ps) {
  if (ps->callbacks->end_multi_pattern_scan) {
    ps->callbacks->end_multi_pattern_scan(ps->invoker_state);
  }
}

void end_ele_cis_element_search(CISML_PARSER_T *ps) {
  if (ps->callbacks->end_cis_element_search) {
    ps->callbacks->end_cis_element_search(ps->invoker_state);
  }
}

/*handle_start_doc {{{*/
void handle_start_doc(void *ctx) {
  CISML_PARSER_T *ps = (CISML_PARSER_T*)ctx;
  DEBUG_MSG(DUMP_VERBOSE, "CISML parser - start of CISML document\n");
  if (ps->callbacks->start_cisml) {
    ps->callbacks->start_cisml(ps->invoker_state);
  }
}
/*}}}*/

/*handle_end_doc {{{*/
void handle_end_doc(void *ctx) {
  CISML_PARSER_T *ps = (CISML_PARSER_T*)ctx;
  DEBUG_MSG(DUMP_VERBOSE, "CISML parser - end of CISML document\n");
  if (ps->callbacks->end_cisml) {
    ps->callbacks->end_cisml(ps->invoker_state);
  }
}
/*}}}*/

/*handle_characters {{{*/
int accept_all(char ch) {
  return 1;
}

int accept_not_isspace(char ch) {
  return !isspace(ch);
}

void handle_characters(void *ctx, const xmlChar *ch, int len) {
  CISML_PARSER_T *ps = (CISML_PARSER_T*)ctx;
  if (ps->state == PS_ERROR) return;
  if (ps->udepth) {//we don't know where we are!
    if (ps->callbacks->characters_unknown) {
      ps->callbacks->characters_unknown(ps->invoker_state, ch, len);
    }
  } else {
    CHARBUF_T *buf = &(ps->characters);
    int (*accept) (char) = ps->characters.accept;
    int i, start, end, accepted;
    if (accept) {
      i = 0;
      while (i < len) {
        //find how many to skip
        for (; i < len; ++i) if (accept((char)ch[i])) break;
        if (i >= len) return; //nothing left
        start = i;
        //find how many to accept
        for (++i; i < len; ++i) if (!(accept((char)ch[i]))) break;
        end = i;
        accepted = end - start;
        //increase buffer size if needed
        if ((buf->pos + accepted) >= buf->size) {
          buf->size = buf->pos + accepted + 1;
          buf->buffer = mm_realloc(buf->buffer, buf->size * sizeof(char));
        }
        //copy over the characters to the buffer
        for (i = start; i < end; ++i) {
          buf->buffer[(buf->pos)++] = (char)ch[i];
        }
        buf->buffer[buf->pos] = '\0';
      }
    }
  }
}
/*}}}*/

/*handle_start_ele {{{*/
#define DO_START_ELE(_expected_,_call_,_transition_,_char_accept_) \
  if (strcmp((char*)name, #_expected_ ) == 0) { \
    ps->state = _transition_; \
    _call_ (ps, attrs); \
    ps->characters.accept = _char_accept_; \
    break; \
  }

#define CHECK_START_ELE(_expected_,_transition_,_char_accept_) \
  if (strcmp((char*)name, #_expected_ ) == 0) { \
    ps->state = _transition_; \
    ps->characters.accept = _char_accept_; \
    break; \
  }

#define IGNORE NULL
#define ALL_CHARS accept_all
#define ALL_BUT_SPACE accept_not_isspace

void handle_start_ele(void *ctx, const xmlChar *name, const xmlChar **attrs) {
  CISML_PARSER_T *ps = (CISML_PARSER_T*)ctx;
  int known, allow_unknown;
  if (ps->state == PS_ERROR) return;
  if (ps->udepth) {//we don't know where we are!
    ps->udepth += 1;
  } else {
    //reset the character buffer to the begining
    ps->characters.pos = 0;
    ps->characters.buffer[0] = '\0';
    known = 1; //assume we can find it
    allow_unknown = 0; //are unknowns allowed in this state?
    switch (ps->state) {
      case PS_START:
        DO_START_ELE(cis-element-search, start_ele_cis_element_search, PS_IN_CIS_ELEMENT_SEARCH, IGNORE);
        known = 0;
        break;
      case PS_IN_CIS_ELEMENT_SEARCH:
        CHECK_START_ELE(program-name, PS_IN_PROGRAM_NAME, ALL_CHARS);
        DO_START_ELE(parameters, start_ele_parameters, PS_IN_PARAMETERS, IGNORE);
        DO_START_ELE(multi-pattern-scan, start_ele_multi_pattern_scan, PS_IN_MULTI_PATTERN_SCAN, IGNORE);
        DO_START_ELE(pattern, start_ele_pattern, PS_IN_PATTERN, IGNORE);
        known = 0;
        break;
      case PS_IN_PARAMETERS:
        allow_unknown = 1;
        CHECK_START_ELE(pattern-file, PS_IN_PATTERN_FILE, ALL_CHARS);
        CHECK_START_ELE(sequence-file, PS_IN_SEQUENCE_FILE, ALL_CHARS);
        CHECK_START_ELE(background-seq-file, PS_IN_BACKGROUND_SEQ_FILE, ALL_CHARS);
        CHECK_START_ELE(pattern-pvalue-cutoff, PS_IN_PATTERN_PVALUE_CUTOFF, ALL_CHARS);
        CHECK_START_ELE(sequence-pvalue-cutoff, PS_IN_SEQUENCE_PVALUE_CUTOFF, ALL_CHARS);
        CHECK_START_ELE(site-pvalue-cutoff, PS_IN_SITE_PVALUE_CUTOFF, ALL_CHARS);
        DO_START_ELE(sequence-filtering, start_ele_sequence_filtering, PS_IN_SEQUENCE_FILTERING, IGNORE);
        known = 0;
        break;
      case PS_IN_MULTI_PATTERN_SCAN:
        allow_unknown = 1;
        DO_START_ELE(pattern, start_ele_pattern, PS_IN_PATTERN, IGNORE);
        known = 0;
        break;
      case PS_IN_PATTERN:
        allow_unknown = 1;
        DO_START_ELE(scanned-sequence, start_ele_scanned_sequence, PS_IN_SCANNED_SEQUENCE, IGNORE);
        known = 0;
        break;
      case PS_IN_SCANNED_SEQUENCE:
        allow_unknown = 1;
        DO_START_ELE(matched-element, start_ele_matched_element, PS_IN_MATCHED_ELEMENT, IGNORE);
        known = 0;
        break;
      case PS_IN_MATCHED_ELEMENT:
        allow_unknown = 1;
        CHECK_START_ELE(sequence, PS_IN_SEQUENCE, ALL_CHARS);
        known = 0;
      default:
        known = 0;
    }
    if (!known) {
      if (allow_unknown) {
        ps->udepth = 1;
        if (ps->callbacks->start_unknown) {
          ps->callbacks->start_unknown(ps->invoker_state, name, attrs);
        }
      } else {
        DEBUG_FMT(HIGH_VERBOSE, "CISML parser encountered illegal tag %s while in state %s\n",(char*)name, state_names[ps->state]);
        ps->state = PS_ERROR;
      }
    }
  }
}
/*}}}*/

/*handle_end_ele {{{*/
#define DO_END_ELE(_expected_,_call_,_state_,_transition_) \
  case _state_: \
    if (strcmp((char*)name, #_expected_) == 0) { \
      known = 1; \
      _call_ (ps); \
      if (ps->state == _state_) ps->state = _transition_; \
    } \
    break
#define CHECK_END_ELE(_expected_,_state_,_transition_) \
  case _state_: \
    if (strcmp((char*)name, #_expected_) == 0) { \
      known = 1; \
      ps->state = _transition_; \
    } \
    break

void handle_end_ele(void *ctx, const xmlChar *name) {
  CISML_PARSER_T *ps = (CISML_PARSER_T*)ctx;
  int known = 0;
  if (ps->state == PS_ERROR) return;
  if (ps->udepth) {
    ps->udepth -= 1; 
    if (ps->callbacks->end_unknown) {
      ps->callbacks->end_unknown(ps->invoker_state, name);
    }
  } else {
    switch (ps->state) {
      DO_END_ELE(cis-element-search, end_ele_cis_element_search, PS_IN_CIS_ELEMENT_SEARCH, PS_END);
      DO_END_ELE(program-name, end_ele_program_name, PS_IN_PROGRAM_NAME, PS_IN_CIS_ELEMENT_SEARCH);
      DO_END_ELE(parameters, end_ele_parameters, PS_IN_PARAMETERS, PS_IN_CIS_ELEMENT_SEARCH);
      DO_END_ELE(pattern-file, end_ele_pattern_file, PS_IN_PATTERN_FILE, PS_IN_PARAMETERS);
      DO_END_ELE(sequence-file, end_ele_sequence_file, PS_IN_SEQUENCE_FILE, PS_IN_PARAMETERS);
      DO_END_ELE(background-seq-file, end_ele_background_seq_file, PS_IN_BACKGROUND_SEQ_FILE, PS_IN_PARAMETERS);
      DO_END_ELE(pattern-pvalue-cutoff, end_ele_pattern_pvalue_cutoff, PS_IN_PATTERN_PVALUE_CUTOFF, PS_IN_PARAMETERS);
      DO_END_ELE(sequence-pvalue-cutoff, end_ele_sequence_pvalue_cutoff, PS_IN_SEQUENCE_PVALUE_CUTOFF, PS_IN_PARAMETERS);
      DO_END_ELE(site-pvalue-cutoff, end_ele_site_pvalue_cutoff, PS_IN_SITE_PVALUE_CUTOFF, PS_IN_PARAMETERS);
      CHECK_END_ELE(sequence-filtering, PS_IN_SEQUENCE_FILTERING, PS_IN_PARAMETERS);
      DO_END_ELE(multi-pattern-scan, end_ele_multi_pattern_scan, PS_IN_MULTI_PATTERN_SCAN, PS_IN_CIS_ELEMENT_SEARCH);
      case PS_IN_PATTERN:
        if (strcmp((char*)name, "pattern") == 0) {
          known = 1;
          if (ps->multi == MP_MULTI_PATTERN) {
            ps->state = PS_IN_MULTI_PATTERN_SCAN;
          } else {
            ps->state = PS_IN_CIS_ELEMENT_SEARCH;
          }
          end_ele_pattern(ps);
        }
        break;
      DO_END_ELE(scanned-sequence, end_ele_scanned_sequence, PS_IN_SCANNED_SEQUENCE, PS_IN_PATTERN);
      DO_END_ELE(matched-element, end_ele_matched_element, PS_IN_MATCHED_ELEMENT, PS_IN_SCANNED_SEQUENCE);
      DO_END_ELE(sequence, end_ele_sequence, PS_IN_SEQUENCE, PS_IN_MATCHED_ELEMENT);
    }
    if (!known) {
      DEBUG_FMT(HIGH_VERBOSE,"CISML parser error: unexpected end tag %s in state %s\n", (char*)name, state_names[ps->state]);
      ps->state = PS_ERROR;
    }
  }
}
/*}}}*/

/*parse_cisml {{{*/
int parse_cisml(CISML_CALLBACKS_T *callbacks, void *state, const char *file_name) {
  CISML_PARSER_T cisml_parser;
  CHARBUF_T *buf;
  xmlSAXHandler handler;
  int result;

  DEBUG_FMT(HIGHER_VERBOSE, "CISML parser processing \"%s\"\n", file_name);

  //initilise parser state
  cisml_parser.callbacks = callbacks;
  cisml_parser.invoker_state = state;
  cisml_parser.state = PS_START;
  cisml_parser.multi = MP_UNDECIDED;
  cisml_parser.udepth = 0;

  //set up character buffer
  buf = &(cisml_parser.characters);
  buf->buffer = mm_malloc(sizeof(char)*10);
  buf->buffer[0] = '\0';
  buf->size = 10;
  buf->pos = 0;

  //set up handler
  memset(&handler, 0, sizeof(xmlSAXHandler));
  handler.startDocument = handle_start_doc;
  handler.endDocument = handle_end_doc;
  handler.characters = handle_characters;
  handler.startElement = handle_start_ele;
  handler.endElement = handle_end_ele;

  //parse
  result = xmlSAXUserParseFile(&handler, &cisml_parser, file_name);

  //clean up memory
  free(buf->buffer);
  buf->buffer = NULL;

  //check result
  if (result != 0) {
    DEBUG_FMT(HIGH_VERBOSE, "CISML parser halted due to SAX error; error code given: %d\n", result);
  } else {
    if (cisml_parser.state == PS_END) {
      DEBUG_MSG(HIGHER_VERBOSE, "CISML parser completed\n");
    } else {
      DEBUG_FMT(HIGH_VERBOSE, "CISML parser did not reach end state; final state was %s\n", state_names[cisml_parser.state]);
    }
  }

  // return true on success
  return (result == 0 && cisml_parser.state == PS_END);
}
/*}}}*/

