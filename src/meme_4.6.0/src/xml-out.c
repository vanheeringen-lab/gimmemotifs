/**************************************************************************
 * FILE: xml-out.c
 * CREATE DATE: 4/May/2010
 * AUTHOR: James Johnson
 * PROJECT: shared
 * COPYRIGHT: UQ, 2010
 * DESCRIPTION: Utility functions for writing XML files.
 **************************************************************************/
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>


/**********************************************************************/
/*
 * scans through a null terminated input string and copies it into the
 * buffer replacing '&' with '&amp;', '<' with '&lt;', '>' with '&gt;' and
 * optionally '"' with '&quot;'. Extra is an optional parameter but if it
 * is specified then it will contain the amount that buffer needs to be
 * expanded by to fit the entire translated input. The buffer is returned
 * and is always null terminated.
 *
 * WARNING: do not use this twice in one printf with the same buffer as the
 * second call will overwrite the buffer before the printf is evaluated.
 */
/**********************************************************************/
char* replace_xml_chars(char *input, char *buffer, int buffer_size, int replace_quote, int *extra) {
  int i, j, k;
  char *copy;
  for (i = 0, j = 0; input[i] != '\0'; ++i) {
    //put the next item to be copied across into the copy string
    switch (input[i]) {
      case '&': copy = "&amp;"; break;
      case '<': copy = "&lt;"; break;
      case '>': copy = "&gt;"; break;
      case '"':
        if (replace_quote) {
          copy = "&quot;"; break;
        }
      //fall through
      default://copy across the character (if possible)
        if (j < buffer_size) buffer[j] = input[i];
        ++j; continue;
    }
    //copy the item into the buffer
    for (k = 0; copy[k] != '\0'; ++k, ++j) {
      if (j < buffer_size) {
        buffer[j] = copy[k];
      }
    }
  }
  if (j >= buffer_size) {
    if (extra) *extra = j - buffer_size + 1;
    if (buffer_size > 0) buffer[buffer_size - 1] = '\0';
  } else {
    if (extra) *extra = 0;
    buffer[j] = '\0';
  }
  return buffer;
} /* replace_xml_chars */


/**********************************************************************/
/*
 * scans through a null terminated input string and copies it into the
 * buffer replacing '&' with '&amp;', '<' with '&lt;', '>' with '&gt;' and
 * optionally '"' with '&quot;'. If the buffer is not large enough it
 * will be expanded. The offset can be used to append to an existing string
 * already in the buffer.
 *
 * WARNING: do not use this twice in one printf with the same buffer as the
 * second call will overwrite the buffer before the printf is evaluated.
 */
/**********************************************************************/
char* replace_xml_chars2(char *input, char **expandable_buffer, int *buffer_size, int offset, int replace_quote) {
  char *expanded_buffer, *out;
  int extra = 0;
  assert(*buffer_size >= 0);
  assert(offset >= 0);
  out = replace_xml_chars(input, (*expandable_buffer)+offset, *buffer_size - offset, replace_quote, &extra);
  if (extra) {
    expanded_buffer = realloc(*expandable_buffer, sizeof(char) * (*buffer_size + extra));
    if (expanded_buffer) {
      *expandable_buffer = expanded_buffer;
      *buffer_size += extra;
    } else {
      //memory allocation failure
      fprintf(stderr, "FATAL: replace_xml_chars2 - realloc failed to expand buffer by %d bytes.\n", extra);
      exit(1);
    }
    out = replace_xml_chars(input, (*expandable_buffer)+offset, *buffer_size - offset, replace_quote, &extra);
  }
  assert(extra == 0);
  return out;
} /* replace_xml_chars2 */
