/******************************************************************************
 * FILE: seq-reader-from-fasta.c
 * AUTHOR: Charles Grant, Bill Noble, Tim Bailey
 * CREATION DATE: 2010-09-17
 * COPYRIGHT: 2010 UW
 * 
 * This file contains the concrete implementation for a 
 * data-block-reader UDT for reading sequence segments from a FASTA file.
 *****************************************************************************/

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alphabet.h"
#include "seq-reader-from-fasta.h"

typedef struct seq_reader_from_fasta {
  BOOLEAN_T at_start_of_line;
  int current_position;
  char* filename;
  size_t filename_len; // Includes trailing '\0'
  size_t filename_buffer_len;
  char* sequence_header;
  size_t sequence_header_len; // Includes trailing '\0'
  size_t sequence_buffer_len;
  FILE *fasta_file;
  char *alphabet;
} SEQ_READER_FROM_FASTA_T;

// Forward declarations

DATA_BLOCK_READER_T *new_seq_reader_from_fasta(const char *filename);
void free_seq_reader_from_fasta(DATA_BLOCK_READER_T *reader);
BOOLEAN_T close_seq_reader_from_fasta(DATA_BLOCK_READER_T *reader);
BOOLEAN_T reset_seq_reader_from_fasta(DATA_BLOCK_READER_T *reader);
BOOLEAN_T read_seq_header_from_seq_reader_from_fasta(
  SEQ_READER_FROM_FASTA_T *fasta_reader
);
BOOLEAN_T get_seq_name_from_seq_reader_from_fasta(
  DATA_BLOCK_READER_T *reader, 
  char **name // OUT
);
BOOLEAN_T go_to_next_sequence_in_seq_reader_from_fasta(
  DATA_BLOCK_READER_T *reader
);
BOOLEAN_T get_next_data_block_from_seq_reader_from_fasta(
  DATA_BLOCK_READER_T *reader, 
  DATA_BLOCK_T *data_block
);

/******************************************************************************
 * This function creates an instance of a data block reader UDT for reading
 * sequence segments from a FASTA file.
 *****************************************************************************/
DATA_BLOCK_READER_T *new_seq_reader_from_fasta(const char *filename) {
  SEQ_READER_FROM_FASTA_T *fasta_reader = mm_malloc(sizeof(SEQ_READER_FROM_FASTA_T) * 1);
  fasta_reader->at_start_of_line = TRUE;
  int filename_len = strlen(filename) + 1;
  fasta_reader->filename = mm_malloc(sizeof(char)* filename_len);
  fasta_reader->filename_len = filename_len;
  strncpy(fasta_reader->filename, filename, filename_len);
  fasta_reader->sequence_header = NULL;
  fasta_reader->sequence_header_len = 0;
  fasta_reader->alphabet = get_alphabet(TRUE /* include ambiguity codes */);
  if (
    open_file(
      filename, 
      "r", 
      TRUE, 
      "FASTA", 
      "sequences", 
      &(fasta_reader->fasta_file
    )) == FALSE) {
    die("Couldn't open the file %s.\n", filename);
  }
  fasta_reader->current_position = 0;

  DATA_BLOCK_READER_T *reader = new_data_block_reader(
    (void *) fasta_reader,
    free_seq_reader_from_fasta,
    close_seq_reader_from_fasta,
    reset_seq_reader_from_fasta,
    get_next_data_block_from_seq_reader_from_fasta,
    go_to_next_sequence_in_seq_reader_from_fasta,
    get_seq_name_from_seq_reader_from_fasta
  );
  return reader;
}

/******************************************************************************
 * This function frees an instance of the sequence FASTA reader UDT.
 *****************************************************************************/
void free_seq_reader_from_fasta(DATA_BLOCK_READER_T *reader) {
  SEQ_READER_FROM_FASTA_T *fasta_reader 
    = (SEQ_READER_FROM_FASTA_T *) get_data_block_reader_data(reader);
  myfree(fasta_reader->filename);
  fasta_reader->filename_len = 0;
  fasta_reader->filename_buffer_len = 0;
  myfree(fasta_reader->sequence_header);
  fasta_reader->sequence_header_len = 0;
  fasta_reader->sequence_buffer_len = 0;
  myfree(fasta_reader);
}

/******************************************************************************
 * This function closes a sequence FASTA reader UDT.
 *****************************************************************************/
BOOLEAN_T close_seq_reader_from_fasta(DATA_BLOCK_READER_T *reader) {
  BOOLEAN_T result = FALSE;
  SEQ_READER_FROM_FASTA_T *fasta_reader 
    = (SEQ_READER_FROM_FASTA_T *) get_data_block_reader_data(reader);
  fasta_reader->current_position = 0;
  if (fasta_reader->fasta_file != NULL) {
    if (fclose(fasta_reader->fasta_file) == EOF) {
      die(
        "Error closing file: %s.\nError message: %s\n", 
        fasta_reader->filename, 
        strerror(errno)
      );
    }
    else {
      result = TRUE;
    }
  }
  return result;
}

/******************************************************************************
 * This function resets a sequence FASTA reader UDT.
 *****************************************************************************/
BOOLEAN_T reset_seq_reader_from_fasta(DATA_BLOCK_READER_T *reader) {
  SEQ_READER_FROM_FASTA_T *fasta_reader 
    = (SEQ_READER_FROM_FASTA_T *) get_data_block_reader_data(reader);
  if (fasta_reader->fasta_file == stdin) {
    die("Unable to rewind when reading sequence from standard input\n");
  }
  else {
    rewind(fasta_reader->fasta_file);
  }
  fasta_reader->current_position = -1;
  fasta_reader->at_start_of_line = TRUE;
  return TRUE;
}

/******************************************************************************
 * This function reads the entire sequence header at the start of a new sequence.
 * The current position is assumed to be start of a new sequence.
 * Read from the current position to the end of the current line.
 *
 * Returns TRUE if it was able to read the sequence text, FALSE if 
 * EOF reached before the terminal newline was found. Dies if other errors
 * are encountered.
 *****************************************************************************/
BOOLEAN_T read_seq_header_from_seq_reader_from_fasta(
  SEQ_READER_FROM_FASTA_T *fasta_reader
) {

  int result = FALSE;

  // Initial allocation of sequence buffer
  const size_t initial_buffer_len = 100;
  if (fasta_reader->sequence_header == NULL) {
     fasta_reader->sequence_header = mm_malloc(sizeof(char) * initial_buffer_len);
     fasta_reader->sequence_buffer_len = initial_buffer_len;
  }

  // Look for EOL
  int c = 0;
  int seq_index = 0;
  while((c = fgetc(fasta_reader->fasta_file)) != EOF) {

    if (seq_index >= fasta_reader->sequence_buffer_len) {
      // Need to grow buffer
      fasta_reader->sequence_header
        = mm_realloc(
            fasta_reader->sequence_header, 
            2 * fasta_reader->sequence_buffer_len
          );
      fasta_reader->sequence_buffer_len = 2 * fasta_reader->sequence_buffer_len;
    }

    if (c == '\n') {
      // Found EOL
      fasta_reader->sequence_header[seq_index] = '\0';
      fasta_reader->sequence_header_len = seq_index + 1;
      fasta_reader->at_start_of_line = TRUE;
      result = TRUE;
      break;
    }
    else {
      // Keep looking for EOL
      fasta_reader->sequence_header[seq_index] = c;
      ++seq_index;
    }

  }

  // At this point c is EOL or EOF
  if (c == EOF) {
    if (ferror(fasta_reader->fasta_file)) {
      // EOF could actually indicate an error.
      die(
        "Error reading file:%s.\nError message: %s\n", 
        fasta_reader->filename,
        strerror(ferror(fasta_reader->fasta_file))
      );
    }
    else if (feof(fasta_reader->fasta_file)) {
        // Reached EOF before reaching EOL for the sequence.
        fasta_reader->sequence_header[0] = '\0';
        fasta_reader->sequence_header_len = 0;
    }
  }

  return result;

}

/******************************************************************************
 * This function gets the name of the current sequence from the data block
 * reader. The name of the sequence is passed using the name parameter.
 * The caller is responsible for freeing the memory for the sequence name.
 *
 * Returns TRUE if successful, FALSE if there is no current sequence, as 
 * at the start of the file.
 *****************************************************************************/
BOOLEAN_T get_seq_name_from_seq_reader_from_fasta(
  DATA_BLOCK_READER_T *reader, 
  char **name // OUT
) {
  BOOLEAN_T result = FALSE;
  SEQ_READER_FROM_FASTA_T *fasta_reader 
    = (SEQ_READER_FROM_FASTA_T *) get_data_block_reader_data(reader);
  if (fasta_reader->sequence_header == NULL || fasta_reader->sequence_header_len <= 0) {
    result = FALSE;
  }
  else {
    int name_len = 0;
    for (name_len = 0; name_len < fasta_reader->sequence_header_len; ++name_len) {
      if (isspace(fasta_reader->sequence_header[name_len])) {
          break;
      }
    }
    myassert(
      TRUE, 
      name_len <= fasta_reader->sequence_header_len, 
      "Error parsing seq.  name\n"
    );
    char *buffer = mm_malloc(sizeof(char) * (name_len + 1));
    strncpy(buffer, fasta_reader->sequence_header, name_len);
    buffer[name_len] = 0;
    *name = buffer;
    result = TRUE;
  }

  return result;
}

/******************************************************************************
 * Read from the current position in the file to the first prior after the
 * start of the next sequence. Set the value of the current sequence.
 *
 * Returns TRUE if it was able to advance to the next sequence, FALSE if 
 * EOF reached before the next sequence was found. Dies if other errors
 * encountered.
 *****************************************************************************/
BOOLEAN_T go_to_next_sequence_in_seq_reader_from_fasta(DATA_BLOCK_READER_T *reader) {
  BOOLEAN_T result = FALSE;
  SEQ_READER_FROM_FASTA_T *fasta_reader 
    = (SEQ_READER_FROM_FASTA_T *) get_data_block_reader_data(reader);
  fasta_reader->current_position = 0;
  int c = 0;
  while((c = fgetc(fasta_reader->fasta_file)) != EOF) {
    if (fasta_reader->at_start_of_line == TRUE && c == '>') {
      break;
    }
    else if (c == '\n') {
      fasta_reader->at_start_of_line = TRUE;
    }
    else {
      fasta_reader->at_start_of_line = FALSE;
    }
  }
  // At this point c is '>' or EOF
  if (c == '>') {
    result = read_seq_header_from_seq_reader_from_fasta(fasta_reader);
  }
  else {
    if (ferror(fasta_reader->fasta_file)) {
      die(
        "Error reading file:%s.\nError message: %s\n", 
        fasta_reader->filename,
        strerror(ferror(fasta_reader->fasta_file))
      );
    }
    else if (feof(fasta_reader->fasta_file)) {
        // Reached EOF before reaching the start of the sequence
        result = FALSE;
    }
  }
  return result;
}

/******************************************************************************
 * Fills in the next data block for the sequence. 
 * During the first call for the sequence it fills in the full data block.
 * On successive calls, shifts the sequence in the block down one position
 * and reads one more character.
 * 
 * Returns TRUE if it was able to completely fill the block, FALSE if 
 * the next sequence or EOF was reached before the block was filled.
 * Dies if other errors encountered.
 *****************************************************************************/
BOOLEAN_T get_next_data_block_from_seq_reader_from_fasta(
  DATA_BLOCK_READER_T *reader, 
  DATA_BLOCK_T *data_block
) {

  BOOLEAN_T result = FALSE;
  char *raw_seq = get_sequence_from_data_block(data_block);
  int block_size = get_block_size_from_data_block(data_block);
  int num_read = get_num_read_into_data_block(data_block);

  SEQ_READER_FROM_FASTA_T *fasta_reader 
    = (SEQ_READER_FROM_FASTA_T *) get_data_block_reader_data(reader);

  if (num_read == block_size) {
    // Block is alread full, shift all elements in the block down by one position
    // FIXME CEG: Inefficient, replace with circular buffer.
    memmove(raw_seq, raw_seq + 1, block_size - 1);
    num_read = block_size - 1;
    raw_seq[num_read] = 0;
  }

  int c = 0;
  while((c = fgetc(fasta_reader->fasta_file)) != EOF) {
    if (isspace(c)) {
      // Skip over white space
      if (c == '\n') {
        fasta_reader->at_start_of_line = TRUE;
      }
      else {
        fasta_reader->at_start_of_line = FALSE;
      }
      continue;
    }
    else if (c == '>' && fasta_reader->at_start_of_line == TRUE) {
      // We found the start of a new sequence while trying
      // to fill the block. Leave the block incomplete.
      c = ungetc(c, fasta_reader->fasta_file);
      if (ferror(fasta_reader->fasta_file)) {
        die(
          "Error while reading file:%s.\nError message: %s\n", 
          fasta_reader->filename,
          strerror(ferror(fasta_reader->fasta_file))
        );
      }
      raw_seq[num_read] = 0;
      break;
    }
    else {
      // Fill in another character in the block
      raw_seq[num_read] = toupper(c);
      // Check that character is legal in alphabet. 
      // If not, replace with wild card character.
      if (!char_in_string(fasta_reader->alphabet, raw_seq[num_read])) {
        ALPH_T alphabet_type = which_alphabet();
        char wild_card = (alphabet_type ==  PROTEIN_ALPH ? ANY_AMINO : ANY_BASE);
        raw_seq[num_read] = wild_card;
	      fprintf(
          stderr, 
          "Warning: %c is not a valid character in alphabet %s.\n"
          "         Converting %c to %c.\n",
		      c,
          fasta_reader->alphabet,
		      c,
          raw_seq[num_read]
        );
      }
      ++num_read;
      if (num_read == block_size) {
        // block is full
        result = TRUE;
        break;
      }
    }
  }

  if (c == EOF && ferror(fasta_reader->fasta_file)) {
    die(
      "Error while reading file:%s.\nError message: %s\n", 
      fasta_reader->filename,
      strerror(ferror(fasta_reader->fasta_file))
    );
  }

  ++fasta_reader->current_position;
  set_start_pos_for_data_block(data_block, fasta_reader->current_position);
  set_num_read_into_data_block(data_block, num_read);
  return result;

}
