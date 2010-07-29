
#include <stdio.h>

//note: matrix.h must be declared before array.h
#include "matrix.h"
#include "alphabet.h"
#include "array.h"
#include "config.h"
#include "motif.h"
#include "ceqlogo.h"
#include "io.h"
#include "meme-io.h"
#include "utils.h"
#include "user.h"


#define DEFAULT_PSEUDOCOUNTS 0
#define LOGO_PREFIX "logo"


VERBOSE_T verbosity = QUIET_VERBOSE;


void copy_and_sanatise_name(char *target, char *source, int max_len) {
  int i;
  i = 0;
  while (i < max_len && source[i] != '\0') {
    if ((source[i] >= 'a' && source[i] <= 'z') || 
        (source[i] >= 'A' && source[i] <= 'Z') || 
        (source[i] >= '0' && source[i] <= '9') ||
        source[i] == '_') {
      target[i] = source[i];
    } else {
      target[i] = '_';
    }
    ++i;
  }
  if (i < max_len) {
    target[i] = '\0';
  }
}

void generate_ceq_logos(char *meme_path, char *output_dir) {
  int i, num_motifs, dir_len, prefix_len, path_len;
  ARRAY_T *background;
  BOOLEAN_T has_reverse_strand;
  char *path, *alphabet;
  double logo_height, logo_width;
  /* this is a strange way to do this, 
   * but meme-io forces me to preallocate 
   * the number of motifs I want, even 
   * when I have no idea how many there are! */
  MOTIF_T motifs[MAX_MOTIFS];  

  logo_height = LOGOHEIGHT;
  //make the path
  dir_len = strlen(output_dir);
  prefix_len = strlen(LOGO_PREFIX);
  path_len = dir_len + 1 + prefix_len + MAX_MOTIF_ID_LENGTH + 1;
  path = malloc(sizeof(char)*path_len);
  strncpy(path, output_dir, path_len);
  if (path[dir_len-1] != '/') {
    path[dir_len] = '/';
    path[++dir_len] = '\0';
  }
  strncpy(path+dir_len, LOGO_PREFIX, path_len - dir_len);

  // Read all motifs into an array.
  read_meme_file(meme_path,
		 NULL, // bg file name
		 DEFAULT_PSEUDOCOUNTS,
     REQUIRE_PSPM,
		 &num_motifs,
		 motifs, 
		 NULL,//motif occurrences
		 &has_reverse_strand,
		 &background);

  // global alphabet is set by read_meme_file
  alphabet = get_alphabet(FALSE);

  if (create_output_directory(output_dir, TRUE, (verbosity >= NORMAL_VERBOSE))) {
    // Failed to create output directory.
    exit(1);
  }

  for(i = 0; i < num_motifs; i++) {
    logo_width = motifs[i].length <= MAXLOGOWIDTH ? motifs[i].length : MAXLOGOWIDTH;
    copy_and_sanatise_name(path+(dir_len+prefix_len), motifs[i].id, path_len - (dir_len + prefix_len)); 
    CL_create2(
      &(motifs[i]), 			// motif
      "", 			          // no title 
      NULL, 			        // no second motif
      "", 			          // no x-axis label
      FALSE, 			        // no error bars
      FALSE,			        // ssc
      logo_height,		    // logo height (cm)
      logo_width,		      // logo width (cm)
      alphabet, 	        // alphabet
      0, 			            // no offset to second motif
      path,			          // output file path
      "MEME (no SSC)"		  // program name
    );
  }
  for(i=0; i<num_motifs; i++) { // free all loaded motifs
    free_motif(&motifs[i]);
  }
  free_array(background); // not used 
  free(path);
}



int main(int argc, char** argv) {
  if (argc != 3) {
    fprintf(stderr, "Expected:\nmeme2images <meme file> <output directory>\n");
    exit(1);
  }
  generate_ceq_logos(argv[1], argv[2]);
  return 0;
}
