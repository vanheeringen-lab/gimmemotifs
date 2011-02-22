/**************************************************************************
 * FILE: xml-out.h
 * CREATE DATE: 4/May/2010
 * AUTHOR: James Johnson
 * PROJECT: shared
 * COPYRIGHT: UQ, 2010
 * DESCRIPTION: Utility functions for writing XML files.
 **************************************************************************/
#ifndef XML_OUT_H
#define XML_OUT_H

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
char* replace_xml_chars(char *input, char *buffer, int buffer_size, int replace_quote, int *extra);


/**********************************************************************/
/*
 * scans through a null terminated input string and copies it into the
 * buffer replacing '&' with '&amp;', '<' with '&lt;', '>' with '&gt;' and
 * optionally '"' with '&quot;'. If the buffer is not large enough it
 * will be expanded.
 *
 * WARNING: do not use this twice in one printf with the same buffer as the
 * second call will overwrite the buffer before the printf is evaluated.
 */
/**********************************************************************/
char* replace_xml_chars2(char *input, char **expandable_buffer, int *buffer_size, int offset, int replace_quote);

#endif
