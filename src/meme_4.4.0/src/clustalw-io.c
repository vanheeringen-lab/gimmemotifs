/****************************************************************************
 * FILE: clustalw-io.c
 * AUTHOR: William Stafford Noble, Charles E. Grant
 * CREATE DATE: 10/27/2004
 * PROJECT: EVOMCAST
 * DESCRIPTION: Read an alignment from a CLUSTAL file (.aln) into memory.
 * COPYRIGHT: 1998-2004, UCSD, UCSC, UW
 ****************************************************************************/
#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include "alignment.h"
#include "clustalw-io.h"
#include "seq.h"
#include "string-list.h"
#include "utils.h"

#define LONGEST_LINE 10000  // Longest input line allowed.

/****************************************************************************
 * Read and parse the header line from a CLUSTALW file.
 * We ignore the remainder of the line.
 *
 * Returns TRUE if the correct header line was found, FALSE otherwise.
 ****************************************************************************/
static BOOLEAN_T read_clustalw_header(FILE* clustalw_file)
{
  BOOLEAN_T result = FALSE;
  static char line[LONGEST_LINE];

  assert(clustalw_file != NULL);

  // Read the CLUSTAL header line.
  if (fgets(line, LONGEST_LINE-1, clustalw_file) != NULL) {

    // The first line of a CLUSTALW file should always
    // start with the words "CLUSTAL " or "CLUSTALW"
    if ((strncmp(line, "CLUSTAL ", 8) == 0) ||
        (strncmp(line, "CLUSTALW", 8) == 0)) {
      result = TRUE;

      /***********************************************************
      // We are not going to parse the remainder of the line
      //
      // Parse the header line. It should consist of four segments:
      // program name, program type, version and a description, 
      // each separated by white space.

      int length = 0;
      char* segments[4];
      int i = 0;
      int segment_index = 1;

      // We've already found the name
      segments[0] = line;
      line[7] = '\0';

      for (segment_index = 1; segment_index < 4; segment_index++) {
        // Find next white space
        while (i < length && !isspace(line[i])) {
          i++;
        }
        if (i >= length) {
          result = FALSE;
          break;
        }
        // Skip over white space
        while (i < length && isspace(line[i])) {
          line[i] = '\0';
          i++;
        }
        if (i >= length) {
          result = FALSE;
          break;
        }
        // We've found the type
        segments[segment_index] = &(line[i]);
      }
      if (result == TRUE) {
        copy_string(name, segments[0]);
        if (name == NULL) result = FALSE;
        copy_string(type, segments[1]);
        if (type == NULL) result = FALSE;
        copy_string(version, segments[2]);
        if (version == NULL) result = FALSE;
        copy_string(description, segments[3]);
        if (description == NULL) result = FALSE;
      }
      **********************************************************/
    }
  }

  return(result);
}

/****************************************************************************
 * Parse a line of sequence data from an alignment block in a CLUSTALW file
 * Adds sequence name and sequence string to the caller provided string lists.
 *
 * Returns TRUE if parse was successful, FALSE otherwise.
 ****************************************************************************/
BOOLEAN_T parse_clustalw_sequence_line(
    char *line,
    int length,
    STRING_LIST_T* names,
    STRING_LIST_T* sequences,
    int*  num_columns
)
{
  BOOLEAN_T result = TRUE;
  char *segments[3];
  int i = 0; 
  int name_index = 0;
  int num_names = 0; 
  int segment_index = 0;
  int sequence_offset = 0;

  assert(line != NULL);
  assert(names != NULL);
  assert(sequences != NULL);
  assert(num_columns != NULL);

  if (isspace(line[0])) {
    // A sequence line can't start with white space
    result = FALSE;
  } else {

    // Parse the sequence line. It should consist of two
    // or three segments: sequence name, sequence string,
    // and an optional residue count. The segments are separated
    // by white space.

    for (segment_index = 0; segment_index < 3; segment_index++) {
      // Skip over white space
      while (i < length && isspace(line[i])) {
        line[i] = '\0';
        i++;
      }
      // Have we gone too far?
      if (i >= length) {
        result = FALSE;
        break;
      }
      if (segment_index == 1) {
        // Record the starting position of the sequence string
        sequence_offset = i;
      }
      segments[segment_index] = &(line[i]);
      // Find next white space
      while (i < length && !isspace(line[i])) {
        i++;
      }
      if (segment_index == 1) {
        // We'll use our own count of columns in preference to
        // the optional field of the clustalw file sequence record
        *num_columns = i - sequence_offset + 1;
      }
      // Have we gone too far?
      if (i >= length) {
        if (segment_index == 0) result = FALSE;
        break;
      }
    }

    if (result == TRUE) {
      // We parsed the line, search for the
      // sequence name in the list of known names
      num_names = get_num_strings(names);
      char *query_name = segments[0];
      int query_name_len = strlen(query_name);
      for (name_index = 0; name_index < num_names; name_index++) {
        char *target_name = get_nth_string(name_index, names);
        int target_name_len = strlen(target_name);
        if (query_name_len == target_name_len 
            && strncmp(query_name, target_name, query_name_len) == 0) {
            break;
        }
      }
      if (name_index == num_names) {
        // Name was not found in list of names 
        // so add name and sequence to lists.
        add_string(segments[0], names);
        add_string(segments[1], sequences);
      } else {
        // This name was already in the list.
        // Append to existing sequences string.
        append_to_nth_string(segments[1], name_index, sequences);
      }
    }
  }
  return result;
}

/****************************************************************************
 * Read and parse an alignment block from a CLUSTALW file
 * Sequence names and strings are stored in the caller provided string lists.
 * Caller must free memory allocated to the consensus string.
 *
 * Returns 1 if a block was read, 0 if EOF was reached,
 * and -1 if an error occured.
 ****************************************************************************/
static int read_clustalw_block(
   FILE* clustalw_file,
   STRING_LIST_T* sequence_names,
   STRING_LIST_T* sequence_strings,
   char**  consensus_string
)
{
  BOOLEAN_T got_block = FALSE;
  BOOLEAN_T in_block = FALSE;
  static char line[LONGEST_LINE]; 
  int num_columns = 0;
  int offset = 0;
  int num_sequences = 0;
  int length = 0;
  int result;

  assert(clustalw_file != NULL);
  assert(sequence_names != NULL);
  assert(sequence_strings != NULL);
  assert(consensus_string != NULL);

  while (fgets(line, LONGEST_LINE - 1, clustalw_file) != NULL) {
    // Remove trailing newline character
    length = strlen(line) - 1;
    if (line[length] == '\n') {
      line[length] = '\0';
    }
    // Skip empty lines
    if (length == 0) {
      continue; 
    }
    if (isspace(line[0])) {
      if (in_block) {
        // A non-empty line with leading whitespace after we've started
        // a sequence block indicates we've reached the consensus string.
        // We need to match the length of the consensus string to the
        // length of the sequence strings.
        offset =  length - num_columns + 1;
        copy_string(consensus_string, line + offset);
        if (*consensus_string != NULL) { 
          // We've successfully reached the end of the block
          got_block = TRUE;
        }
      } else {
        // This file is not in the correct format
      }
      break;
    } else {
      // Non-whitespace should indicate a line of sequence data.
      // Each alignment block line will have two or three segments:
      // a sequence name, a sequence string, and optionally a column count.
      if (parse_clustalw_sequence_line(line, length, sequence_names, 
            sequence_strings, &num_columns)) {
        // The line parsed
        in_block = TRUE;
        num_sequences++;
      } else {
        // The line didn't parse
        got_block = FALSE;
        break;
      }
    }
  }
  // Did we get the block?
  if (got_block) {
    result = 1;
  } else {
    // Was there an error or was it just EOF?
    if (feof(clustalw_file)) {
      result = 0;
    } else {
      // There must have been an error
      result = -1;
    }
  }

  return result;
}

/****************************************************************************
 * Read an alignment from a file in CLUSTALW (.aln) format.
 * If call is successful caller is responsible for freeing alignment object.
 *
 * Return: Was an alignment successfully read?
 ****************************************************************************/
BOOLEAN_T read_clustalw(
   FILE*     clustalw_file,
   ALIGNMENT_T** alignment
)
{
  BOOLEAN_T result = TRUE;
  char* name = NULL;
  char* description = NULL;
  char *consensus_fragment = NULL;
  int i;
  int got_block;
  STRING_LIST_T* names = NULL;
  STRING_LIST_T* sequences = NULL;
  STRING_LIST_T* consensus = NULL;
  SEQ_T** aligned_sequences = NULL;

  assert(clustalw_file != NULL);
  assert(alignment != NULL);

  result = read_clustalw_header(clustalw_file);   
  if (result) {
    names = new_string_list();
    sequences = new_string_list();
    consensus = new_string_list();
    add_string("", consensus);
    // Read until there are no more blocks
    while ((got_block = read_clustalw_block(clustalw_file, names, sequences, 
          &consensus_fragment)) == 1) {
      append_to_nth_string(consensus_fragment, 0, consensus);
      myfree(consensus_fragment);
      consensus_fragment = NULL;
    }
    if (got_block == -1) {
      // An error occured while trying to read a clustal block.
      fprintf(stderr, "Error reading CLUSTAL W block.\n");
      result = FALSE;
    } else {
      // Build the array of sequences for the alignment.
      aligned_sequences = (SEQ_T**) mm_malloc(get_num_strings(sequences) 
          * sizeof(SEQ_T *));
      for (i = 0; i < get_num_strings(names); i++) {
        aligned_sequences[i] = allocate_seq(get_nth_string(i, names), "", 0,
            get_nth_string(i, sequences));
      }
      *alignment = allocate_alignment(name, description, 
        i, aligned_sequences, get_nth_string(0, consensus));
      // Having created the alignment we can free
      // the intermediate components.
      for (i = 0; i < get_num_strings(names); i++) {
        free_seq(aligned_sequences[i]);
      }
    }
    // Clean up intermediate objects as needed
    if (consensus != NULL) free_string_list(consensus);
    if (sequences != NULL) free_string_list(sequences);
    if (names != NULL) free_string_list(names);
    if (aligned_sequences != NULL) myfree(aligned_sequences);
  } else {
    fprintf(stderr, "Error reading CLUSTALW header.\n");
  }

  return result;
}

ALIGNMENT_T* read_alignment_from_clustalw_file(char* filename) {

  BOOLEAN_T result = FALSE;
  char* align_name = NULL;
  FILE* align_file = NULL;
  ALIGNMENT_T* alignment = NULL;

  result = open_file(
    filename, 
    "r", 
    TRUE, 
    "alignment", 
    "alignment", 
    &align_file
  );
  if (result == FALSE) {
    die("Couldn't open the file %s.\n", filename);
  }

  // File should be a CLUSTALW alignment file.
  if (read_clustalw(align_file, &alignment) == TRUE) {
    align_name = strrchr(filename, '/');
    if (align_name) {
      align_name++;
    } else {
      // No '/' in name
      align_name = filename;
    }
  } else {
    die("Error while reading %s as a CLUSTALW file.\n", filename);
  }
  set_alignment_name(align_name, alignment);
  fclose(align_file);

  return alignment;
}

/****************************************************************************
 * Print an alignment in CLUSTALW format.
 ****************************************************************************/
void print_clustalw(
   FILE* outfile,
   BOOLEAN_T show_residue_count,
   ALIGNMENT_T* an_alignment
)
{
  int i = 0;
  int length = 0;
  int num_sequences = 0;
  int seq_index = 0;
  int num_blocks = 0;
  int block_index = 0;
  int remainder = 0;
  int name_length = 0;
  int num_spaces = 0;
  int max_seq_name = 0;
  SEQ_T* seq = NULL;
  char* name = NULL;
  char* consensus_string = NULL;
  char* seq_string = NULL;
  int* residue_counts = NULL;

  char seq_block[BLOCKSIZE + 1];

  assert(an_alignment != NULL);

  // Allocate space to keep track of residue counts
  // for each sequence
  num_sequences = get_num_aligned_sequences(an_alignment);
  residue_counts = (int *) mm_malloc(num_sequences * sizeof(int));
  if (residue_counts == NULL) {
    die("Error allocating memory to track residue counts\n");
  }
  memset(residue_counts, 0, num_sequences * sizeof(int));

  // Calculate how many blocks we must
  // divide the sequences into
  length = get_alignment_length(an_alignment);
  num_blocks = length/BLOCKSIZE;
  remainder = (length % BLOCKSIZE);
  if (remainder) num_blocks++;
  
  // Print the heading line.
  fputs("CLUSTAL W\n", outfile);
  fputs("\n\n", outfile);

  // Calcuate how long the longest sequence name is
  for (seq_index = 0; seq_index < num_sequences; seq_index++) {
    seq = get_alignment_sequence(seq_index, an_alignment);
    name_length = strlen(get_seq_name(seq));
    max_seq_name = max_seq_name > name_length ? max_seq_name : name_length;
  }
  // Add 12 more spaces for to separate the columns
  max_seq_name += 12;

  for (block_index = 0; block_index < num_blocks; block_index++) {
    for (seq_index = 0; seq_index < num_sequences; seq_index++) {
      seq = get_alignment_sequence(seq_index, an_alignment);
      assert(seq != NULL);
      name = get_seq_name(seq);
      fputs(name, outfile);
      // Pad name with white space so that columns align
      num_spaces = max_seq_name - strlen(name);
      for (i = 0; i < num_spaces; i++) {
        fputc(' ', outfile);
      }
      seq_string = get_raw_sequence(seq);
      seq_block[0] = '\0';
      strncat(seq_block, (seq_string + (block_index * BLOCKSIZE)), BLOCKSIZE);
      residue_counts[seq_index] += count_residues(seq_block);
      fputs(seq_block, outfile);
      if (show_residue_count) {
        fprintf(outfile, " %d", residue_counts[seq_index]);
      }
      fputc('\n', outfile);
    }
    seq_block[0] = '\0';
    consensus_string = get_consensus_string(an_alignment);
    strncat(seq_block, consensus_string + (block_index * BLOCKSIZE), BLOCKSIZE);
    for (i = 0; i < max_seq_name; i++) {
      fputc(' ', outfile);
    }
    fprintf(outfile, "%s\n\n", seq_block);
  }
  if (residue_counts) myfree(residue_counts);
}

#ifdef MAIN
#include "simple-getopt.h"

VERBOSE_T verbosity = NORMAL_VERBOSE;

/*************************************************************************
 * int main
 *************************************************************************/
int main(int argc, char *argv[]) {
  char*     clustalw_filename = NULL;
  FILE*     clustalw_file;

  clustalw_filename = argv[1];

  // Open the file for reading.
  if (open_file(clustalw_filename, "r", TRUE, "CLUSTAL", "alignment", 
        &clustalw_file) == 0) {
    exit(1);
  }

  // Read the alignment
  ALIGNMENT_T* alignment = NULL;
  BOOLEAN_T result = read_clustalw(clustalw_file, &alignment);
  if (result) {
    print_clustalw(stdout, TRUE, alignment);
    print_clustalw(stdout, FALSE, alignment);
    free_alignment(alignment);
  } else {
    fprintf(stderr, "Unable to parse an alignment from the file %s.\n",
        clustalw_filename);
  }

  fclose(clustalw_file);
  return(!result);
}
#endif
