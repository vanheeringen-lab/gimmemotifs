/*
 * motif_regexp.c
 *
 *  Created on: 22/09/2008
 *      Author: rob
 */

#include "motif_regexp.h"

/********************************************************************
  The nucleic acid codes supported are:

        A --> adenosine           M --> A C (amino)
        C --> cytidine            S --> G C (strong)
        G --> guanine             W --> A T (weak)
        T --> thymidine           B --> G T C
        U --> uridine             D --> G A T
        R --> G A (purine)        H --> A C T
        Y --> T C (pyrimidine)    V --> G C A
        K --> G T (keto)          N --> A G C T (any)

  The accepted amino acid codes are:

    A  alanine                         P  proline
    B  aspartate or asparagine         Q  glutamine
    C  cystine                         R  arginine
    D  aspartate                       S  serine
    E  glutamate                       T  threonine
    F  phenylalanine                   U  any
    G  glycine                         V  valine
    H  histidine                       W  tryptophan
    I  isoleucine                      Y  tyrosine
    K  lysine                          Z  glutamate or glutamine
    L  leucine                         X  any
    M  methionine
    N  asparagine						(From alphabet.h)
 ********************************************************************/

#define MAX_MOTIF_WIDTH 100

void read_regexp_file(
   char*      filename,          // Name of MEME file  IN
   int*       num_motifs,             // Number of motifs retrieved  OUT
   MOTIF_T*   motifs                 // The retrieved motifs - NOT ALLOCATED!
) {
	FILE*      motif_file;         // MEME file containing the motifs.
	char motif_name[MAX_MOTIF_ID_LENGTH+1];
	char motif_regexp[MAX_MOTIF_WIDTH];
	ARRAY_T* these_freqs;
	MOTIF_T* m;
	int i;

	//Set things to the defaults.
	*num_motifs = 0;

	// Open the given MEME file.
	if (open_file(filename, "r", TRUE, "motif", "motifs", &motif_file) == 0)
		exit(1);

	//Set alphabet - ONLY supports dna.
	set_alphabet(verbosity, "ACGT");

	while (fscanf(motif_file, "%s\t%s", motif_name, motif_regexp) == 2) {
		/*
		 * Now we:
		 * 1. Fill in new motif (preallocated)
		 * 2. Assign name
		 * 3. Convert regexp into frequency table.
		 */

		m = &(motifs[*num_motifs]);
		set_motif_id(motif_name, m);
		m->length = strlen(motif_regexp);
		/* Store the alphabet size in the motif. */
		m->alph_size = get_alph_size(ALPH_SIZE);
		m->ambigs = get_alph_size(AMBIG_SIZE);
		/* Allocate memory for the matrix. */
		m->freqs = allocate_matrix(m->length, get_alph_size(ALL_SIZE));

		//Set motif frequencies here.
		for (i=0;i<strlen(motif_regexp);i++) {
			switch(toupper(motif_regexp[i])) {
			case 'A':
				set_matrix_cell(i,alphabet_index('A',get_alphabet(TRUE)),1,m->freqs);
				break;
			case 'C':
				set_matrix_cell(i,alphabet_index('C',get_alphabet(TRUE)),1,m->freqs);
				break;
			case 'G':
				set_matrix_cell(i,alphabet_index('G',get_alphabet(TRUE)),1,m->freqs);
				break;
			case 'T':
				set_matrix_cell(i,alphabet_index('T',get_alphabet(TRUE)),1,m->freqs);
				break;
			case 'U':
				set_matrix_cell(i,alphabet_index('U',get_alphabet(TRUE)),1,m->freqs);
				break;
			case 'R': //purines
				set_matrix_cell(i,alphabet_index('G',get_alphabet(TRUE)),1,m->freqs);
				set_matrix_cell(i,alphabet_index('A',get_alphabet(TRUE)),1,m->freqs);
				break;
			case 'Y': //pyramidines
				set_matrix_cell(i,alphabet_index('T',get_alphabet(TRUE)),1,m->freqs);
				set_matrix_cell(i,alphabet_index('C',get_alphabet(TRUE)),1,m->freqs);
				break;
			case 'K': //keto
				set_matrix_cell(i,alphabet_index('G',get_alphabet(TRUE)),1,m->freqs);
				set_matrix_cell(i,alphabet_index('T',get_alphabet(TRUE)),1,m->freqs);
				break;
			case 'M': //amino
				set_matrix_cell(i,alphabet_index('A',get_alphabet(TRUE)),1,m->freqs);
				set_matrix_cell(i,alphabet_index('C',get_alphabet(TRUE)),1,m->freqs);
				break;
			case 'S': //strong
				set_matrix_cell(i,alphabet_index('G',get_alphabet(TRUE)),1,m->freqs);
				set_matrix_cell(i,alphabet_index('C',get_alphabet(TRUE)),1,m->freqs);
				break;
			case 'W': //weak
				set_matrix_cell(i,alphabet_index('A',get_alphabet(TRUE)),1,m->freqs);
				set_matrix_cell(i,alphabet_index('T',get_alphabet(TRUE)),1,m->freqs);
				break;
			case 'B':
				set_matrix_cell(i,alphabet_index('G',get_alphabet(TRUE)),1,m->freqs);
				set_matrix_cell(i,alphabet_index('T',get_alphabet(TRUE)),1,m->freqs);
				set_matrix_cell(i,alphabet_index('C',get_alphabet(TRUE)),1,m->freqs);
				break;
			case 'D':
				set_matrix_cell(i,alphabet_index('G',get_alphabet(TRUE)),1,m->freqs);
				set_matrix_cell(i,alphabet_index('A',get_alphabet(TRUE)),1,m->freqs);
				set_matrix_cell(i,alphabet_index('T',get_alphabet(TRUE)),1,m->freqs);
				break;
			case 'H':
				set_matrix_cell(i,alphabet_index('A',get_alphabet(TRUE)),1,m->freqs);
				set_matrix_cell(i,alphabet_index('C',get_alphabet(TRUE)),1,m->freqs);
				set_matrix_cell(i,alphabet_index('T',get_alphabet(TRUE)),1,m->freqs);
				break;
			case 'V':
				set_matrix_cell(i,alphabet_index('G',get_alphabet(TRUE)),1,m->freqs);
				set_matrix_cell(i,alphabet_index('C',get_alphabet(TRUE)),1,m->freqs);
				set_matrix_cell(i,alphabet_index('A',get_alphabet(TRUE)),1,m->freqs);
				break;
			case 'N':
				set_matrix_cell(i,alphabet_index('A',get_alphabet(TRUE)),1,m->freqs);
				set_matrix_cell(i,alphabet_index('C',get_alphabet(TRUE)),1,m->freqs);
				set_matrix_cell(i,alphabet_index('G',get_alphabet(TRUE)),1,m->freqs);
				set_matrix_cell(i,alphabet_index('T',get_alphabet(TRUE)),1,m->freqs);
				break;
			}
		}

	    /* Compute values for ambiguous characters. */
		for (i = 0; i < m->length; i++) {
		    these_freqs = get_matrix_row(i, m->freqs);
		    fill_in_ambiguous_chars(FALSE, these_freqs);
		}

		/* Compute and store the motif complexity. */
		m->complexity = compute_motif_complexity(m);

		//Move our pointer along to do the next motif.
		(*num_motifs)++;
	}
}

