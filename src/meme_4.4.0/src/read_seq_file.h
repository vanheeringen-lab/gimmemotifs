/*
 * $Id: read_seq_file.h 776 2006-05-10 17:27:13Z cegrant $
 * 
 * $Log$
 * Revision 1.1  2005/07/29 19:09:14  nadya
 * Initial revision
 *
 */

#ifndef READ_DATA_FILE_H
#define READ_DATA_FILE_H

extern DATASET *read_seq_file(
  char *file_name,		/* name of file to open */
  char *alpha,			/* alphabet used in sequences */
  BOOLEAN use_comp,		/* use complementary strands, too */
  double seqfrac 		/* fraction of input sequences to use */
);
extern SAMPLE *get_sample_by_name(
  char *sample_name
);

#endif
