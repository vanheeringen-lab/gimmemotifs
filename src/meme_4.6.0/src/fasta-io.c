/****************************************************************************
 * FILE: fasta-io.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 9/14/98
 * PROJECT: MHMM
 * DESCRIPTION: Read sequences from a FASTA file into memory.
 * COPYRIGHT: 1998-2008, UCSD, UCSC, UW
 ****************************************************************************/
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include "utils.h"
#include "alphabet.h"
#include "io.h"
#include "seq.h"
#include "fasta-io.h"

/****************************************************************************
 * Find the beginning of the next sequence, and read the sequence name
 * and description. The name is everything after the '>' to the
 * first white space. The description is everything after all white
 * space following the name to the end of the line.
 * After calling *name will point to a buffer containing the name, 
 * and *description will point to a buffer containing the description. 
 * The caller is * responsible for freeing the memory in these buffers.
 *
 * The return value is TRUE if a title line was read and FALSE otherwise.
 ****************************************************************************/
static BOOLEAN_T read_title_line
  (FILE* fasta_file,
   char** name,
   char** description)
{
  char *id_line = NULL;  // Line containing the ID and comment.

  // Read until the first occurrence of ">".
  do {
    id_line = getline2(fasta_file);
  } while (id_line != NULL && id_line[0] != '>');

  if (id_line == NULL) {
    // We never found a title line
    return FALSE;
  }

  // id_line will contain trailing '\n'
  size_t line_length = strlen(id_line) - 1;
  id_line[line_length] = 0;
  if (line_length == 1) {
    // The sequence must have a name.
    // but the line contained only a '>'
    die("Error reading sequence ID.\n");
  }

  // Separate sequence name from description
  size_t name_start = 1;
  size_t description_start = line_length;
  size_t name_length = 0;
  size_t description_length = 0;
  size_t line_index = 0;

  // Find end of sequence name (ends at first white space)
  for (line_index = 1; line_index < line_length; ++line_index) {
    if (isspace(id_line[line_index])) {
      description_start = line_index;
      break;
    }
    name_length = line_index;
  }

  // Skip over white space at the start of description
  for(; description_start < line_length; ++description_start) {
    if (!isspace(id_line[description_start])) {
      break;
    }
  }

  if (description_start < line_length) {
    description_length = line_length - description_start;
  }

  *name = mm_malloc(sizeof(char) * (name_length + 1));
  strncpy(*name, id_line + 1, name_length);
  (*name)[name_length] = 0;
  if (description_length > 0) {
    *description = mm_malloc(sizeof(char) * (description_length + 1));
    strncpy(*description, id_line + description_start, description_length);
    (*description)[description_length] = 0;
  }
  else {
    *description = mm_malloc(sizeof(char));
    (*description)[0] = '\0';
  }

  myfree(id_line);
  return(TRUE);
}

/****************************************************************************
 * Read raw sequence until a '>' is encountered or too many letters
 * are read.  The new sequence is appended to the end of the given
 * sequence.
 *
 * Return: Was the sequence read completely?
 ****************************************************************************/
static BOOLEAN_T read_raw_sequence
  (FILE* fasta_file,   // Input Fasta file.
   char* name,         // Sequence ID (used in error messages).
   int   max_chars,    // Maximum number of allowed characters.
   char* raw_sequence) // Pre-allocated sequence.
{
  // char a_char;
  // tlb; change a_char to integer so it will compile on SGI
  int a_char;
  int i_seq;
  BOOLEAN_T return_value = TRUE;

  // Start at the end of the given sequence.
  i_seq = strlen(raw_sequence);
  assert((int)strlen(raw_sequence) < max_chars);

  // Get the basic alphabet
  char *alphabet = get_alphabet(FALSE);

  // Set the character to replace unknowns with
  char any_char = ANY_AMINO;
  if (strlen(alphabet) == NUM_BASES) any_char = ANY_BASE;

  // Get the full alphabet (with ambigs)
  alphabet = get_alphabet(TRUE);

  // Read character by character.
  while ((a_char = getc(fasta_file)) != EOF) {

    // Check for the beginning of the next sequence.
    if (a_char == '>') {
      // Put the ">" back onto the stream for the next call to find.
      ungetc(a_char, fasta_file);
      break;
    }

    // Skip non-alphabetic characters.
    if (!isalpha((int)a_char)) {
      if ((a_char != ' ') && (a_char != '\t') && (a_char != '\n')) {
	fprintf(stderr, "Warning: Skipping character %c in sequence %s.\n",
		a_char, name);
      }

    } else {

      // Convert invalid characters to any_char.
      a_char = toupper((int)a_char);
      if (!char_in_string(alphabet, a_char)) {
	fprintf(stderr, "Warning: Converting illegal character %c to %c ",
		a_char, any_char);
	fprintf(stderr, "in sequence %s.\n", name);
	a_char = any_char;
      }
      raw_sequence[i_seq] = a_char;
      i_seq++;
    }
    if (i_seq >= max_chars) {
      return_value = FALSE;
      break;
    }
  }
  raw_sequence[i_seq] = '\0';

  return(return_value);
}

/****************************************************************************
 * Read one sequence from a file in Fasta format.
 *
 * Return: Was a sequence successfully read?
 ****************************************************************************/
BOOLEAN_T read_one_fasta
  (FILE*     fasta_file,
   int       max_seq_length,
   SEQ_T**   sequence)
{
  char *name = NULL;     // Just the sequence ID.
  char *desc = NULL;     // Just the comment field.
  char* buffer = NULL;         // The sequence, as it's read in.

  // Read the title line.
  if (!read_title_line(fasta_file, &name, &desc)) {
    return(FALSE);
  }
  
  // Create local static storage space.
  if (buffer == NULL) {
    buffer = mm_malloc(sizeof(char) * max_seq_length);
  }

  // Read the sequence.
  buffer[0] = '\0';
  if (!read_raw_sequence(fasta_file, name, max_seq_length, buffer)) {
    die("Sequence %s is too long (maximum=%d).\n", name, max_seq_length);
  }
    
  // Allocate the new sequence object.
  *sequence = allocate_seq(name, desc, 0, buffer);
  myfree(name);
  myfree(desc);
  myfree(buffer);

  // Tell the user what's up.
  if (verbosity >= DUMP_VERBOSE) {
    fprintf(stderr, 
	    "Read sequence %s of length %d.\n",
	    get_seq_name(*sequence), get_seq_length(*sequence));
  }
  return(TRUE);
}

/****************************************************************************
 * Read up to N letters of one sequence from a file in Fasta format.
 * If the given sequence is NULL, start from the beginning of the next
 * sequence.  Otherwise, append to the given sequence.
 *
 * Return: Was a sequence segment successfully read?
 ****************************************************************************/
BOOLEAN_T read_one_fasta_segment
  (FILE*     fasta_file,
   int       max_chars,
   SEQ_T**   sequence)
{
  char *name;     // Just the sequence ID.
  char *desc;     // Just the comment field.
  int         offset = 0;             // Offset of this sequence.
  char*       new_sequence = (char*)mm_calloc(max_chars + 1, sizeof(char));
  BOOLEAN_T   is_complete;  // Is this a complete sequence?

  // Extract information from the current sequence.
  if (*sequence != NULL) {
    // FIXME CEG This is inefficient. We should really add a 
    // method to allow setting the sequence in SEQ_T objects
    char *old_name = get_seq_name(*sequence);
    int len = strlen(old_name) + 1;
    name = mm_malloc(sizeof(char) * len);
    strcpy(name, old_name);
    char *old_desc = get_seq_description(*sequence);
    len = strlen(old_desc) + 1;
    desc = mm_malloc(sizeof(char) * len);
    strcpy(desc, old_desc);
    assert(get_seq_length(*sequence) < max_chars);
    strcpy(new_sequence, get_raw_sequence(*sequence));
    offset = get_seq_offset(*sequence);
    free_seq(*sequence);
    *sequence = NULL;
  }
  
  // Read information from the file.
  else {
    if (!read_title_line(fasta_file, &name, &desc)) {
      myfree(new_sequence);
      return(FALSE);
    }
  }

  // Read until we run out of space.
  is_complete = read_raw_sequence(fasta_file, name, max_chars,
				  new_sequence);

  // Store the new sequence.
  *sequence = allocate_seq(name, desc, offset, new_sequence);
  set_complete(is_complete, *sequence);
  myfree(name);
  myfree(desc);
  myfree(new_sequence);

  return(TRUE);
}

/****************************************************************************
 * Read all the sequences from a FASTA file at once.
 ****************************************************************************/
#define NUM_ALLOC 1000 /* Allocate how many sequences at once? */
void read_many_fastas
  (FILE*      fasta_file,
   int        max_seq, // Maximum sequence length.
   int*       num_seqs,
   SEQ_T***   sequences)
{
  int i_seq;         /* Index of the current sequence. */
  int num_allocated; /* Number of pointers currently allocated. */

  /* Allocate initial memory. */
  num_allocated = NUM_ALLOC;
  *sequences = (SEQ_T**)mm_malloc(sizeof(SEQ_T*) * num_allocated);

  /* Read the sequences one by one. */
  i_seq = 0;
  while (read_one_fasta(fasta_file, max_seq, &((*sequences)[i_seq]))) {
    i_seq++;

    /* Allocate more space, if need be. */
    if (i_seq >= num_allocated) {
      num_allocated += NUM_ALLOC;
      *sequences = (SEQ_T**)mm_realloc(*sequences, 
				      sizeof(SEQ_T*) * num_allocated);
    }
  }

  /* Record the total number of sequences. */
  *num_seqs = i_seq;

  /* Complain if nothing was read. */
  if (i_seq == 0) {
    die("Failed to read a single sequence from the given FASTA file.\n");
  }
}

/****************************************************************************
 * Print a single sequence in FASTA format.
 ****************************************************************************/
#define FASTA_LINE 50
void print_fasta
  (SEQ_T*    sequence,
   FILE*     outfile)  
{
  int   i_seq;
  int   seq_length;
  char* raw_sequence;

  fprintf(outfile, ">%s %s\n", get_seq_name(sequence),
	  get_seq_description(sequence));

  raw_sequence = get_raw_sequence(sequence);
  i_seq = 0;
  seq_length = strlen(raw_sequence);
  while (seq_length - i_seq > FASTA_LINE) {
    fprintf(outfile, "%.*s\n", FASTA_LINE, &(raw_sequence[i_seq]));
    i_seq += FASTA_LINE;
  }
  fprintf(outfile, "%s\n\n", &(raw_sequence[i_seq]));
}

/****************************************************************************
 *  Create a sequence object from the first sequence in a FASTA file.
 ****************************************************************************/
SEQ_T* read_sequence_from_file(char* filename) {
  
  BOOLEAN_T result = FALSE;
  FILE* seq_file = NULL;
  SEQ_T* seq = NULL;

  set_alphabet(verbosity, "ACGT");

  result = open_file(
    filename, 
    "r", 
    1, 
    "sequence", 
    "sequence", 
    &seq_file
  );
  if (result != TRUE) {
    die("Couldn't open the file %s.\n", filename);
  }
  if (read_one_fasta(seq_file, MAX_SEQ, &seq) == FALSE) {
    die("Error while reading %s as a FASTA file.\n", filename);
  } 
  if (seq_file != NULL) {
    fclose(seq_file);
  }

  return seq;
}

#ifdef MAIN
#include "simple-getopt.h"

VERBOSE_T verbosity = NORMAL_VERBOSE;

/*************************************************************************
 * int main
 *************************************************************************/
int main(int argc, char *argv[])
{
  char*     fasta_filename = NULL;
  FILE*     fasta_file;
  SEQ_T*    sequence = NULL;
  SEQ_T**   sequences;
  int       num_seqs;
  int       i_seq;
  BOOLEAN_T many = FALSE;
  BOOLEAN_T dna = FALSE;
  int       max_chars = -1;

  // Define command line options.
  cmdoption const options[] = {
    {"many", NO_VALUE},
    {"dna", NO_VALUE},
    {"blocksize", REQUIRED_VALUE},
    {"verbosity", REQUIRED_VALUE}
  };
	int option_count = 4;
  int option_index = 0;

  // Define the usage message.
  char      usage[1000] = "";
  strcat(usage, "USAGE: fasta-io [options] <sequences>\n");
  strcat(usage, "\n");
  strcat(usage, "   Options:\n");
  strcat(usage, "     -many\n");
  strcat(usage, "     -dna\n");
  strcat(usage, "     -blocksize <n>\n");
  strcat(usage, "     -verbosity 1|2|3|4|5 (default=2)\n");
  strcat(usage, "\n");

	simple_setopt(argc, argv, option_count, options);

  // Parse the command line.
  while (1) { 
    int c = 0;
    char* option_name = NULL;
    char* option_value = NULL;
		const char *  message = NULL;

    // Read the next option, and break if we're done.
    c = simple_getopt(&option_name, &option_value, &option_index);
    if (c == 0) {
      break;
    } else if (c < 0) {
    	simple_getopterror(&message);
      die("Error processing command line options (%s)\n", message);
    }

    if (strcmp(option_name, "many") == 0) {
      many = TRUE;
    } else if (strcmp(option_name, "dna") == 0) {
      dna = TRUE;
    } else if (strcmp(option_name, "blocksize") == 0) {
      max_chars = atoi(option_value);
    } else if (strcmp(option_name, "verbosity") == 0) {
      verbosity = (VERBOSE_T)atoi(option_value);
    }
  }

  // Read the single required argument.
  if (option_index + 1 != argc) {
    fprintf(stderr, "%s", usage);
    exit(1);
  }
  fasta_filename = argv[option_index];

  if (dna) {
    set_alphabet(verbosity, "ACGT");
  } else {
    set_alphabet(verbosity, "ACDEFGHIKLMNPQRSTVWY");
  }

  // Open the file for reading.
  if (open_file(fasta_filename, "r", TRUE, "FASTA", "sequences", &fasta_file)
      == 0) {
    exit(1);
  }

  if (many) {

    // Read all the sequences into memory.
    read_many_fastas(fasta_file, MAX_SEQ, &num_seqs, &sequences);

    // Print each one and then free it.
    for (i_seq = 0; i_seq < num_seqs; i_seq++) {
      print_fasta(sequences[i_seq], stdout);
      free_seq(sequences[i_seq]);
    }
    myfree(sequences);


  } else if (max_chars == -1) {

    // Read each sequence one at a time, print and free it.
    while (read_one_fasta(fasta_file, MAX_SEQ, &sequence)) {
      print_fasta(sequence, stdout);
      free_seq(sequence);
    }
  } else {

    // Read each sequence one block at a time, print and free it.
    while (read_one_fasta_segment(fasta_file, max_chars, &sequence)) {
      print_fasta(sequence, stdout);
      if (is_complete(sequence)) {
	free_seq(sequence);
	sequence = NULL;
      } else {
	shift_sequence(get_seq_length(sequence), sequence);
      }
    }
  }
    
  fclose(fasta_file);
     
  return(0);
}
#endif


