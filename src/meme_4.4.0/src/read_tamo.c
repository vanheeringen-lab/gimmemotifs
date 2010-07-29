
/********************************************************************
 * FILE: read_tamo.c
 * AUTHOR: Fabian Buske
 * CREATE DATE: 07/07/2008
 * PROJECT: MEME suite
 * COPYRIGHT: 2008, UQ
 *
 * 	Read a tamo motif file.
 *
 * 	Supports the tamo format.
 * 	http://jura.wi.mit.edu/fraenkel/regcode/release_v24/index.html#MotifFormat
 ********************************************************************/

#include <assert.h>
#include <ctype.h>
#include "read_tamo.h"
#include "string-list.h"
#include "alphabet.h"
#include "matrix.h"
#include <stdlib.h>
#include <stdio.h>
#include "meme-io.h"

// Nucleotide alphabest order as in motif.h
extern char alphabet[];

/**********************************************************************
  read_tamo

  Reads in a tamo motif file
**********************************************************************/
BOOLEAN read_tamo(
  char*		tamo_file,			// file containing motifs IN
  double 	pseudocount,		// pseudocount (in frequency)
  int*      num_motifs,         // Number of motifs retrieved  OUT
  MOTIF_T** motif               // The retrieved motifs OUT
)
{
  assert(tamo_file != NULL);
  int i,j,k;
  char c = ' ';
  char* entry = NULL;

  // open tamo file
  FILE *data_file;
  data_file = fopen(tamo_file, "r");  	/* Open in TEXT mode */

  set_alphabet(verbosity, "ACGT");
  int alph_size = get_alph_size(ALPH_SIZE);

  MOTIF_T* motifs = NULL;
  BOOLEAN_T readSequences = FALSE;
  char* motif_source = NULL;
  STRING_LIST_T* lines = new_string_list();

  for (i=0,(*num_motifs)=0; (c=fgetc(data_file))!=EOF; ) {
	if (c == '#'){
		i=0;
		Skip_eol(c, data_file);		/* jump to end of line */
	} else if (c== '\n' || c== '\r') {
    	Resize(entry, i+1, char);
    	entry[i] = '\0';
    	if (readSequences){
    		// create motif
			if (entry[0]=='*'){ // process motif
				readSequences = FALSE;
				Skip_eol(c, data_file);	/* jump to end of line */

				// check if all sequences have same length
	    		BOOLEAN_T check = TRUE;
	    		int len = 0;
	    		for (j=0;j<get_num_strings(lines);++j){
	    			if (len != 0 && strlen(get_nth_string(j,lines)) != len){
	    				fprintf(stderr,"Could not read in motif %s. Skip this instance",motif_source);
	    				check = FALSE;
	    				break;
	    			}
	    			len = strlen(get_nth_string(j,lines));
	    		}

	    		if (check){
	    			Resize(motifs, (*num_motifs)+1, MOTIF_T);
					MATRIX_T* matrix = allocate_matrix(len,get_alph_size(ALL_SIZE));
					for (k=0;k<len;++k){
						double a_freq=0.0;
						double c_freq=0.0;
						double g_freq=0.0;
						double t_freq=0.0;
						for (j=0;j<get_num_strings(lines);++j){
							if (toupper(get_nth_string(j,lines)[k])=='A')
								++a_freq;
							else if (toupper(get_nth_string(j,lines)[k])=='C')
								++c_freq;
							else if (toupper(get_nth_string(j,lines)[k])=='G')
								++g_freq;
							else if (toupper(get_nth_string(j,lines)[k])=='T')
								++t_freq;
						}
						ARRAY_T* counts = get_matrix_row(k, matrix);
						int index = 0;
						index = alphabet_index('A', alphabet);
						set_array_item(index, (a_freq+pseudocount)/get_num_strings(lines), counts);
						index = alphabet_index('C', alphabet);
						set_array_item(index, (c_freq+pseudocount)/get_num_strings(lines), counts);
						index = alphabet_index('G', alphabet);
						set_array_item(index, (g_freq+pseudocount)/get_num_strings(lines), counts);
						index = alphabet_index('T', alphabet);
						set_array_item(index, (t_freq+pseudocount)/get_num_strings(lines), counts);

						// Normalize the first alph_size positions. (MEME prints six
						// digits of precision).
						normalize_subarray(0, alph_size, 0.00001, counts);
						/* Compute values for ambiguous characters. */
						fill_in_ambiguous_chars(FALSE, counts);
					} // motif frequences
					MOTIF_T* m = allocate_motif(motif_source,matrix);
					motifs[*num_motifs] = *m;
					free_matrix(matrix);

				    // Compute and store the motif complexity.
				    (motifs[*num_motifs]).complexity = compute_motif_complexity(&(motifs[*num_motifs]));
				    ++(*num_motifs);

	    		}
				myfree(motif_source);
				free_string_list(lines);
				lines = new_string_list();
			} // end process motif
			else { // read motif instance
				char* instance = NULL;
				int l=0;
    			while (l < strlen(entry)){
    				if (entry[l]==' ')
    					break;
    				else{
	    				if ((l % TCHUNK) == 0) {
	    					Resize(instance, l+TCHUNK, char);
	    				}
	    				instance[l] = entry[l];
	    				++l;
    				}
    			}
    			Resize(instance, l+1, char);
    			instance[l] = '\0';
				add_string(instance, lines);
				myfree(instance);
			} // end motif instance
    	} // motif sequences
    	else { // look for next motif
	    	// check if motif
	    	char* term = "Motif";
	    	BOOLEAN_T match = FALSE;
	    	if (strlen(entry)>= strlen(term)){
	    		match = TRUE;
	    		for (i=0;i<strlen(term);++i){
	    			if (entry[i] != term[i]){
	    				match = FALSE;
	    				break;
	    			}
	    		}
	    	}
	    	if (match){
	    		readSequences = TRUE;
	    	} else { /* check if source */
	    		term = "Source:";
	    		int l;
	    		if (strlen(entry)>= strlen(term)){
		    		match = TRUE;
		    		for (l=0;l<strlen(term);++l){
		    			if (entry[l] != term[l]){
		    				match = FALSE;
		    				break;
		    			}
		    		}
		    	}
	    		if (match){
	    			l = strlen(term);
	    			int count = 0;
	    			while (l < strlen(entry)){
	    				if ((count % TCHUNK) == 0) {
	    					Resize(motif_source, count+TCHUNK, char);
	    				}
	    				motif_source[count++] = entry[l];
	    				++l;
	    			}
	    			Resize(motif_source, count+1, char);
	    			motif_source[count] = '\0';
	    		}
	    	}
    	} // end look for next motif
    	myfree(entry);
    	i=0;
	} else {
		if ((i % TCHUNK) == 0) {
			Resize(entry, i+TCHUNK, char);
		}
		entry[i++] = c;			/* non-blank: add to name */
	}
  }

  if (lines != NULL)
	  free_string_list(lines);

  *motif = motifs;
  fclose(data_file);
  return TRUE;
} /* read_tamo */
