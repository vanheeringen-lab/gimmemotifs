
#include <libxml/parser.h>

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

#include "binary-search.h"
#include "red-black-tree.h"
#include "utils.h"

/*debugging macros {{{*/
VERBOSE_T verbosity = QUIET_VERBOSE;
BOOLEAN_T status = TRUE;
time_t status_last = 0; //last time a status message was output
time_t status_delay = 5; //minimum time between status messages

#define DISPLAY_STATUS(status_msg_format, ...) { \
  if (status) { \
    fprintf(stderr, status_msg_format, __VA_ARGS__); \
    status_last = time(NULL); \
  } \
} \

#define SKIPABLE_STATUS(status_msg_format, ...) { \
  if (status) { \
    time_t status_new = time(NULL); \
    if (status_new > (status_last + status_delay)) { \
      fprintf(stderr, status_msg_format, __VA_ARGS__); \
      status_last = status_new; \
    } \
  } \
} 

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

/* state constants - PS stands for Parser State {{{*/
#define PS_ERROR 0
#define PS_START 1
#define PS_IN_MAST 2
#define PS_IN_MODEL 21
#define PS_IN_COMMAND_LINE 2101
#define PS_IN_MAX_CORRELATION 2102
#define PS_IN_REMOVE_CORRELATED 2103
#define PS_IN_STRAND_HANDLING 2104
#define PS_IN_TRANSLATE_DNA 2105
#define PS_IN_MAX_SEQ_EVALUE 2106
#define PS_IN_ADJ_HIT_PVALUE 2107
#define PS_IN_MAX_HIT_PVALUE 2108
#define PS_IN_MAX_WEAK_PVALUE 2109
#define PS_IN_HOST 2110
#define PS_IN_WHEN 2111
#define PS_IN_ALPHABET 22
#define PS_IN_LETTER 221
#define PS_IN_MOTIFS 23
#define PS_IN_MOTIF 231
#define PS_IN_CORRELATION 232
#define PS_IN_NOS 233
#define PS_IN_EXPECT 2331
#define PS_IN_SEQUENCES 24
#define PS_IN_DATABASE 241
#define PS_IN_SEQUENCE 242
#define PS_IN_SCORE 2421
#define PS_IN_SEG 2422
#define PS_IN_DATA 24221
#define PS_IN_HIT 24222
#define PS_IN_RUNTIME 25
#define PS_END 3
/*}}}*/

/* structs and enums {{{*/
typedef enum {Combine, Separate, Norc, Protein} STYPE;
typedef enum {PRESET, BGFILE, SEQCOMP} ATYPE;

typedef struct alphsym ALPHSYM_T;
struct alphsym {
  char symbol;
  int ambig;
  double bg_value;
};

typedef struct alphabet ALPHABET_T;
struct alphabet {
  int dna;
  int source;
  char *bg_file;
  int count;
  ALPHSYM_T *symbols;
};

typedef struct database DATABASE_T;
struct database {
  char *id;
  int num;
  char *source;
  char *name;
  char *last_mod_date;
  int seq_count;
  long residue_count;
  int is_dna;
};

typedef struct databases DATABASES_T;
struct databases {
  int count;
  DATABASE_T **all_databases; //in file order
  RBTREE_T *lookup; //by id
  int total_hit_count;
};

typedef struct motif MOTIF_T;
struct motif {
  char *id;
  int num;
  char *name;
  int width;
  int ws; //width in seq
  char *best_f;
  char *best_r;
  int bad; //has the motif been excluded due to a bad correlation?
  double *correlations; //correlations with those of num less than our num, array of length num, correlation with motif num x is at position x-1.
};

typedef struct motifs MOTIFS_T;
struct motifs {
  char *source;
  char *name;
  char *last_mod_date;
  int count;
  MOTIF_T **all_motifs; //in file order
  RBTREE_T *lookup; //by id
  char *diagram; //norminal order and spacing diagram
  int diagram_motifs;
  int diagram_remain;
};

typedef struct score SCORE_T;
struct score {
  double combined_pvalue;
  double evalue;
  int db_num;
  int seq_num;
  int strand;
  int frame;
  long file_loc;
};

typedef struct hit HIT_T;
struct hit {
  long pos;
  long gap;
  int strand;
  MOTIF_T *motif;
  double pvalue;
  char *match;
  char *translated;
};

typedef struct seg SEG_T;
struct seg {
  long start;
  char *data;
  int length;
  HIT_T **hits;
  int count;
};

typedef struct sequence SEQUENCE_T;
struct sequence {
  char *id;
  DATABASE_T *db;
  int num;
  char *name;
  char *comment;
  long length;
  int known;
  int has_score1;
  SCORE_T score1; //for the forward strand or both
  int has_score2;
  SCORE_T score2; //for the reverse strand
  SEG_T *segs;   //the segs
  int seg_count; //the number of segs allocated
  int seg_pos;   //the next seg position to use
  HIT_T **hits;  //the hits
  int hit_count; //the number of hits allocated
  int hit_pos;   //the next hit position to use
};

typedef struct charbuf CHARBUF_T;
struct charbuf {
  int (*accept) (char);
  char *buffer;
  int pos;
  int size;
};

typedef struct parser_state PARSER_STATE_T;
struct parser_state {
  char *txt_file_name;
  int state;
  int udepth; //how deep are we from a tag we know?
  FILE *save; //temporary store of hits that we can't print yet
  FILE *hit_file; //table of hits
  FILE *diag_file; //table of motif diagrams
  FILE *note_file; //table of annotated sequences
  char *version; //mast version
  char *release; //mast release
  char *command_line;
  double max_correlation;
  int remove_correlated;
  int translate_dna; //will the dna databases need translating to protein for the motifs
  int db_dna; //are the databases dna
  int strand_handling; 
  double max_seq_evalue;
  int adj_hit_pvalue;
  double max_hit_pvalue;
  double max_weak_pvalue;
  char *host;
  char *when;
  double runtime;
  ALPHABET_T alphabet;
  DATABASES_T databases;
  MOTIFS_T motifs;
  SEQUENCE_T current_seq;
  CHARBUF_T characters;
  RBTREE_T *postponed_scores;
};

typedef struct multi MULTI_T;
struct multi {
  int count;
  char **options;
  int *outputs;
  int *target;
};
/*}}}*/

/* output mast text {{{*/

/* constants for printing control */
#define MMSN 34
#define PAGEWIDTH 80
#define MAXID (80 - MMSN - 7 - PVALUE - LENGTH)
#define PVALUE 8
#define LENGTH 6
#define MAXIDLINES 10

/* macro to copy one open file to another; noop if file NULL */
#define copy_file(in_file, out_file) \
  if (in_file != NULL) { \
    int c; \
    rewind(in_file); \
    while ((c = getc(in_file)) != EOF) putc(c, out_file); \
  }

/*
 * output_mast_txt
 * this takes all that has been loaded and written into separate files
 * and combines it all into the mast text output
 */
void output_mast_txt(PARSER_STATE_T * ps) {
  FILE *mast_out;
  int i, j, bad_motif_count;
  BOOLEAN_T doc = TRUE;
  char *stars = 
"********************************************************************************";
  char *mtype = (ps->translate_dna || !ps->db_dna) ? "peptide" : "nucleotide";
  char *dbtype = (!ps->db_dna) ? "peptide" : "nucleotide";

  mast_out = fopen(ps->txt_file_name, "w");
  if (!mast_out) {
    die("Unable to open \"%s\" for output, error was given as: %s\n", strerror(errno));
  }
  /* 
    announce the program 
  */
  i = strlen(ps->release);
  fprintf(mast_out,"%s\n", stars);
  fprintf(mast_out,"MAST - Motif Alignment and Search Tool\n");
  fprintf(mast_out,"%s\n", stars);
  fprintf(mast_out,
"\tMAST version %s (Release date: %.*s)\n\n"
"\tFor further information on how to interpret these results or to get\n"
"\ta copy of the MAST software please access http://meme.nbcr.net.\n",
    ps->version, i, ps->release);
  fprintf(mast_out,"%s\n", stars);

  /* 
    print reference citation 
  */
  fprintf(mast_out,"\n\n%s\n", stars);
  fprintf(mast_out,"REFERENCE\n");
  fprintf(mast_out,"%s\n", stars);
  fprintf(mast_out,
"\tIf you use this program in your research, please cite:\n"
"\n"
"\tTimothy L. Bailey and Michael Gribskov,\n"
"\t\"Combining evidence using p-values: application to sequence homology\n"
"\tsearches\", Bioinformatics, 14(48-54), 1998.\n"
  );
  fprintf(mast_out,"%s\n", stars);

  /* 
    print info on the database 
  */
  fprintf(mast_out,"\n\n%s\n", stars);
  fprintf(mast_out,"DATABASE AND MOTIFS\n");
  fprintf(mast_out,"%s\n", stars);
  for (i = 0; i < ps->databases.count; ++i) {
    DATABASE_T *db = ps->databases.all_databases[i];
    fprintf(mast_out,"\tDATABASE %s (%s)\n", db->name, dbtype);
    fprintf(mast_out,"\tLast updated on %s\n", db->last_mod_date);
    /* print number of sequences in database */
    fprintf(mast_out,"\tDatabase contains %d sequences, %ld residues\n\n", db->seq_count, db->residue_count);
  }
  /* print handling of DNA strands */
  if (ps->strand_handling == Combine) {
    fprintf(mast_out,"\tScores for positive and reverse complement strands are combined.\n\n");
  } else if (ps->strand_handling == Separate) {
    fprintf(mast_out,"\tPositive and reverse complement strands are scored separately.\n\n");
  } else if (ps->strand_handling == Norc) {
    fprintf(mast_out,"\tReverse complement strands are not scored.\n\n");
  }


  /* 
    print info on the motifs 
  */
  fprintf(mast_out,"\tMOTIFS %s (%s)\n", ps->motifs.name, mtype);
  fprintf(mast_out,"\tMOTIF WIDTH BEST POSSIBLE MATCH\n");
  fprintf(mast_out,"\t----- ----- -------------------\n");
  bad_motif_count = 0;
  for (i = 0; i < ps->motifs.count; i++) {
    MOTIF_T *motif = ps->motifs.all_motifs[i];
    if (motif->bad) {
      bad_motif_count++;
      if (ps->remove_correlated) continue;
    }
    fprintf(mast_out,"\t%3d   %3d   %*.*s\n", motif->num, motif->width, 
      motif->width, motif->width, motif->best_f);
  }
  if (ps->motifs.diagram != NULL) 
    fprintf(mast_out,"\n\tNominal ordering and spacing of motifs:\n\t %s\n", ps->motifs.diagram);

  /* 
   print motif correlations 
  */
  if ((!(ps->remove_correlated) && ps->motifs.count > 1) || (ps->motifs.count - bad_motif_count > 1)) {
    fprintf(mast_out,"\n\tPAIRWISE MOTIF CORRELATIONS:\n");
    fprintf(mast_out,"\tMOTIF");
    for (i = 0; i < ps->motifs.count - 1; i++) {
      MOTIF_T *motif = ps->motifs.all_motifs[i];
      if (ps->remove_correlated && motif->bad) continue;
      fprintf(mast_out," %5d", motif->num);
    }
    fprintf(mast_out,"\n");
    fprintf(mast_out,"\t-----");
    for (i = 0; i < ps->motifs.count - 1; i++) {
      MOTIF_T *motif = ps->motifs.all_motifs[i];
      if (ps->remove_correlated && motif->bad) continue;
      fprintf(mast_out," -----");
    }
    fprintf(mast_out,"\n");
    for (i = 1; i < ps->motifs.count; i++) {                 /* from motif */
      MOTIF_T *from = ps->motifs.all_motifs[i];
      if (ps->remove_correlated && from->bad) continue;
      fprintf(mast_out,"\t  %2d ", from->num);
      for (j=0; j<i; j++) {                     /* to motif */
        MOTIF_T * to = ps->motifs.all_motifs[j];
        if (ps->remove_correlated && to->bad) continue;
        fprintf(mast_out," %5.2f", from->correlations[j]);
      } /* to motif */
      fprintf(mast_out,"\n");
    } /* from motif */
    if (bad_motif_count) {
      fprintf(mast_out,
"\tCorrelations above %4.2f may cause some combined p-values and\n"
"\tE-values to be underestimates.\n", ps->max_correlation 
      );
      if (ps->remove_correlated) {
        fprintf(mast_out,"\tRemoved motif");
      } else {
        fprintf(mast_out,"\tRemoving motif");
      }
      if (bad_motif_count > 1) fprintf(mast_out,"s "); else fprintf(mast_out," ");
      for (i = 0, j = 0; i < ps->motifs.count; i++) {
        MOTIF_T *motif = ps->motifs.all_motifs[i];
        if (motif->bad) {
          if (++j == bad_motif_count) {
            if (j > 1) fprintf(mast_out," and ");
          } else {
            if (j > 1) fprintf(mast_out,", ");
          }
          fprintf(mast_out,"%d", motif->num);
        }
      }
      if (ps->remove_correlated) {
        fprintf(mast_out,
          " because they have correlation > %4.2f\n\twith the remaining motifs.\n", 
          ps->max_correlation);
      } else {
        fprintf(mast_out," from the query may be advisable.\n");
      }
    } else {
      fprintf(mast_out,"\tNo overly similar pairs (correlation > %4.2f) found.\n",
        ps->max_correlation);
    }
  } /* nmotifs > 1 */

  /* print background model frequencies */
  if (ps->alphabet.source != SEQCOMP) {
    int i, pcol;
    char *c;
    fprintf(mast_out,"\n\tRandom model letter frequencies (from %s):",
      ps->alphabet.source == BGFILE ? ps->alphabet.bg_file : "non-redundant database");
    for (i=0, pcol=80; i < ps->alphabet.count; i++) {
      ALPHSYM_T *alphsym = ps->alphabet.symbols+i;
      pcol += 8;          			/* start of printed thing */
      if (pcol >= 80) {pcol=15; fprintf(mast_out,"\n\t");}
      fprintf(mast_out,"%c %5.3f ", alphsym->symbol, alphsym->bg_value);
    }
    fprintf(mast_out,"\n");
  } else {
    fprintf(mast_out,"\n\tUsing random model based on each target sequence composition.\n");
  }
 
  /* end database and motif section */
  fprintf(mast_out,"%s\n", stars);

  /* 
    print table of hits documentation
  */
  fprintf(mast_out,"\n\n%s\nSECTION I: HIGH-SCORING SEQUENCES\n%s\n", stars, stars);
  if (doc) {
    char *fs1 = ps->translate_dna && ps->strand_handling==Separate ? 
      "\n\t- The strand and frame of the (best) motif match(es) is shown." : 
      ps->translate_dna ? "\n\t- The frame of the (best) motif match(es) is shown." : "";
    char *fs2 = ps->translate_dna ? 
      "\n\t  Frames 1, 2, and 3 are labeled a, b c, respectively." : "";
    fprintf(mast_out,
"\t- Each of the following %d sequences has E-value less than %g.\n",
      ps->databases.total_hit_count, ps->max_seq_evalue);
    /* //rank has been permantly disabled
    if (rank > 1) {
      fprintf(mast_out,
"\t- The %d best-matching sequences have been omitted.\n", rank-1);
    }
    */
    fprintf(mast_out,
"\t- The E-value of a sequence is the expected number of sequences\n"
"\t  in a random database of the same size that would match the motifs as\n"
"\t  well as the sequence does and is equal to the combined p-value of the\n"
"\t  sequence times the number of sequences in the database.\n"
"\t- The combined p-value of a sequence measures the strength of the\n"
"\t  match of the sequence to all the motifs and is calculated by\n"
"\t    o finding the score of the single best match of each motif\n"
"\t      to the sequence (best matches may overlap),\n");
    fprintf(mast_out,
"\t    o calculating the sequence p-value of each score,\n"
"\t    o forming the product of the p-values,\n");
    if (ps->motifs.diagram_motifs > 1) {
      fprintf(mast_out,
"\t    o multiplying by the p-value of the observed spacing of\n"
"\t      pairs of adjacent motifs (given the nominal spacing),\n");
    }
    fprintf(mast_out,
"\t    o taking the p-value of the product.\n");
    fprintf(mast_out,
"\t- The sequence p-value of a score is defined as the\n"
"\t  probability of a random sequence of the same length containing\n"
"\t  some match with as good or better a score.\n");
    fprintf(mast_out,
"\t- The score for the match of a position in a sequence to a motif\n"
"\t  is computed by by summing the appropriate entry from each column of\n"
"\t  the position-dependent scoring matrix that represents the motif.\n"
"\t- Sequences shorter than one or more of the motifs are skipped.%s%s\n"
"\t- The table is sorted by increasing E-value.\n", fs1, fs2);
    fprintf(mast_out,"%s\n\n", stars);
  } /* doc */

  /* 
    print table of hits headings and data
  */
  {
    int il = ps->translate_dna || ps->strand_handling==Separate ? MAXID-3 : MAXID+4;/* length of id */
    char *st1 = ps->translate_dna ? "FRAME  " : ps->strand_handling==Separate ? "STRAND " : "";
    char *st2 = ps->translate_dna ? "-----  " : ps->strand_handling==Separate ? "------ " : "";
    char *f = "%-*s %-*s%s %8s %6s\n";
    fprintf(mast_out,f, MMSN, "SEQUENCE NAME", il, "DESCRIPTION", st1,"E-VALUE ","LENGTH");
    fprintf(mast_out,f, MMSN, "-------------", il, "-----------", st2,"--------","------");
    copy_file(ps->hit_file, mast_out);
    fprintf(mast_out,"\n%s\n\n", stars);
  }

  /* 
    print table of diagrams documentation
  */
  fprintf(mast_out,"\n\n%s\nSECTION II: MOTIF DIAGRAMS\n%s\n", stars, stars);
  if (doc) {
    char *fs0 = ps->strand_handling==Separate ? "the strand and\n\t  " : "";
    char *fs1 = ps->strand_handling==Combine ? "s" : "";
    char *fs2 = ps->translate_dna ? "f" : "";			
    char *fs3 = ps->translate_dna ? 
      "\n\t\t    in frame f.  Frames 1, 2, and 3 are labeled a, b c." : ".";
    char *fs4 = ps->strand_handling==Combine ? "\t\t    A minus sign indicates that the occurrence is on the\n\t\t    reverse complement strand.\n" : "";
    char *ptype = ps->adj_hit_pvalue ? "SEQUENCE" : "POSITION";
    char *set = ps->adj_hit_pvalue ?
      "some random subsequence in a set of n,\n where n is the sequence length minus the motif width plus 1," :
      "a single random subsequence of the length of the motif";
    fprintf(mast_out,
"\t- The ordering and spacing of all non-overlapping motif occurrences\n"
"\t  are shown for each high-scoring sequence listed in Section I.\n"
"\t- A motif occurrence is defined as a position in the sequence whose\n"
"\t  match to the motif has %s p-value less than %g.\n"
"\t- The %s p-value of a match is the probability of\n"
"\t  %s\n"
"\t  scoring at least as well as the observed match.\n"
"\t- For each sequence, all motif occurrences are shown unless there\n"
"\t  are overlaps.  In that case, a motif occurrence is shown only if its\n"
"\t  p-value is less than the product of the p-values of the other\n"
"\t  (lower-numbered) motif occurrences that it overlaps.\n"
"\t- The table also shows %sthe E-value of each sequence.\n"
"\t- Spacers and motif occurences are indicated by\n"
"\t   o -d-    `d' residues separate the end of the preceding motif \n"
"\t\t    occurrence and the start of the following motif occurrence\n",
      ptype, ps->max_weak_pvalue, ptype, set, fs0);
    fprintf(mast_out,
"\t   o [%sn%s]  occurrence of motif `n' with p-value less than %g%s\n%s",
      fs1, fs2, ps->max_hit_pvalue, fs3, fs4);
    if (ps->max_weak_pvalue != ps->max_hit_pvalue) fprintf(mast_out,
"\t   o <%sn%s>  occurrence of motif `n' with %g < p-value < %g%s\n%s",
      fs1, fs2, ps->max_hit_pvalue, ps->max_weak_pvalue, fs3, fs4);
    fprintf(mast_out,"%s\n\n", stars);
  } /* doc */

  /*
    print table of diagrams headings and data
  */
  {
    char *st1 = ps->strand_handling==Separate ? "STRAND" : "";	/* strand heading */
    char *st2 = ps->strand_handling==Separate ? "------" : "";	/* strand underline */
    int nl = (*st1=='\0') ? MMSN : MMSN-4;		/* length of name */
    char *f = "%-*s%s %8s  %s\n";			/* format */
    fprintf(mast_out,f, nl, "SEQUENCE NAME", st1, "E-VALUE ", "MOTIF DIAGRAM");
    fprintf(mast_out,f, nl, "-------------", st2, "--------", "-------------");
    copy_file(ps->diag_file, mast_out);
    fprintf(mast_out,"\n%s\n\n", stars);
  }

  /* 
   print table of annotated sequences documentation and data
  */
  if (ps->note_file) {
    fprintf(mast_out,"\n\n%s\n", stars);
    fprintf(mast_out,"SECTION III: ANNOTATED SEQUENCES\n");
    fprintf(mast_out,"%s\n", stars);
    if (doc) {
      char *fs1 = ps->translate_dna ? "/frame" : "";		/* frame string 3 */
      char *fs2 = ps->strand_handling==Combine ? 
          " (a minus sign indicates that\n\t  the occurrence is on the reverse complement strand)" : 
          "";
      char *fs3 = ps->translate_dna && ps->strand_handling==Combine ? 
          " (or its reverse)" : 
        ps->strand_handling==Combine ? 
          " (or its reverse complement)" : 
          "";
      char *fs4 = ps->translate_dna ? "," : ", and";		/* frame string 4 */
      char *fs5 = ps->translate_dna && ps->strand_handling==Combine ? 
	  ", and\n\t   o the protein translation of the match (or its reverse).\n" :
          ps->translate_dna ? 
	  ", and\n\t   o the protein translation of the match.\n" :
          ".\n";
      char *ptype = ps->adj_hit_pvalue ? "SEQUENCE" : "POSITION";
      fprintf(mast_out,
"\t- The positions and p-values of the non-overlapping motif occurrences\n"
"\t  are shown above the actual sequence for each of the high-scoring\n"
"\t  sequences from Section I.\n"
"\t- A motif occurrence is defined as a position in the sequence whose\n"
"\t  match to the motif has %s p-value less than %g as \n"
"\t  defined in Section II.\n"
"\t- For each sequence, the first line specifies the name of the sequence.\n"
"\t- The second (and possibly more) lines give a description of the \n"
"\t  sequence.\n"
"\t- Following the description line(s) is a line giving the length, \n"
"\t  combined p-value, and E-value of the sequence as defined in Section I.\n"
"\t- The next line reproduces the motif diagram from Section II.\n"
"\t- The entire sequence is printed on the following lines.\n"
"\t- Motif occurrences are indicated directly above their positions in the\n"
"\t  sequence on lines showing\n"
"\t   o the motif number%s of the occurrence%s,\n"
"\t   o the position p-value of the occurrence,\n"
"\t   o the best possible match to the motif%s%s\n"
"\t   o columns whose match to the motif has a positive score (indicated \n"
"\t     by a plus sign)%s",
      ptype, ps->max_weak_pvalue, fs1, fs2, fs3, fs4, fs5);
      fprintf(mast_out,"%s\n", stars);
    } /* doc */
    copy_file(ps->note_file, mast_out);
    fprintf(mast_out,"\n%s\n\n", stars);
  } /* annotation table */

  /* display elapsed time */
  fflush(mast_out);
  fprintf(mast_out, "\nCPU: %s\n", ps->host);
  fprintf(mast_out,"Time %f secs.\n\n", ps->runtime);
  /* display the command line */
  fprintf(mast_out, "%s\n", ps->command_line);
}

void print_diagram(char *dia, char *hdr, FILE *file) {
  int j;
  int dia_len = strlen(dia);			/* length of diagram */
  int hlen = strlen(hdr);			/* length of header */

  for (j=0; j < dia_len; ) {
    int remain = dia_len - j;			/* left to print */
    int dlen;					/* room on line */
    char *h = ((j==0) ? hdr : " ");		/* current header */
    dlen = PAGEWIDTH - hlen - 6;
    if (remain <= PAGEWIDTH - hlen) dlen = remain;
    fprintf(file, "%-*.*s%.*s", hlen, hlen, h, dlen, dia+j);
    j += dlen;
    /* continue printing until a good breaking point */
    while (j < dia_len && dia[j-1] != '_' && dia[j-1] != ',') putc(dia[j++], file);
    putc('\n', file);
  }
} /* print_diagram */

/**********************************************************************/
/*
  make_block

  Create a block string:

    [smf(p)] or <smf(p)>
      s   strand (optional)
      m   motif
      f   frame (optional)
      (p)  p-value (optional)  

    returns the number of chars printed 
    or that would be printed if there was
    space (not including null byte). 
    Note that the string is always null 
    terminated (unless size is zero).
*/
/**********************************************************************/
int make_block(
  int m,            /* motif number */
  char *strand,     /* strand */
  int f,            /* frame number; f=0 not translating DNA */
  double thresh,    /* strong motif threshold */
  double p,         /* p-value */
  BOOLEAN_T print_p,/* print p-value in block if TRUE */
  char *block,      /* put block string here */
  int size          /* space left for block string */
)
{
  int printed;
  char left = p < thresh ? '[' : '<';
  char right = p < thresh ? ']' : '>';
  BOOLEAN xlate_dna = (f != 0);
  char *fnames = "abc";          /* frame 1=a, 2=b, 3=c*/

  if (print_p) {          /* print p-value */
    char *bfmt = f ? "%c%s%d%c(%8.2e)%c" : "%c%s%d(%8.2e)%c";
    if (xlate_dna) {          /* str., motif, frame */
      printed = snprintf(block, size, bfmt, left, strand, m, fnames[f-1], p, right);
    } else {            /* strand, motif */
      printed = snprintf(block, size, bfmt, left, strand, m, p, right);
    }
  } else {            /* don't print p-value */
    char *bfmt = f ? "%c%s%d%c%c" : "%c%s%d%c";
    if (xlate_dna) {          /* str., motif, frame */
      printed = snprintf(block, size, bfmt, left, strand, m, fnames[f-1], right);
    } else {            /* strand, motif */
      printed = snprintf(block, size, bfmt, left, strand, m, right);
    }
  }
  return printed;
} /* make_block */


#define DIAG_TEMP_SIZE 100
/*
 * make_diagram
 * generates a motif hit block diagram
 */
char* make_diagram(int strand_handling, int translate_dna, double max_hit_pvalue, BOOLEAN_T print_p, long length, int strand, HIT_T **hits, int count) {
  HIT_T *last_hit;
  char tmp[DIAG_TEMP_SIZE];
  char *diagram;
  int i, mode, printed, pos, remaining, size, strand_hits, first_gap;

  //calculate the first spacer and how many hits we have
  first_gap = length;
  if (strand == 0) {
    strand_hits = count;
    if (count > 0) first_gap = hits[0]->gap;
  } else {
    strand_hits = 0;
    for (i = 0; i < count; ++i) {
      if (hits[i]->strand != strand) continue;
      if (strand_hits++ == 0) first_gap = hits[i]->gap;
    }
  }
  i = 0;
  mode = 1;
  size = 0;
  diagram = NULL;
  remaining = DIAG_TEMP_SIZE;
  if (first_gap != 0) {
    printed = snprintf(tmp, remaining, "%ld", (count == 0 ? length : first_gap));
    if (printed >= remaining) die("make_diagram: tmp buffer too small\n");
  } else {
    printed = 0;
  }
  pos = printed;
  remaining -= printed;
  last_hit = NULL;
  do {
    if (i == count && last_hit != NULL) {
      HIT_T *h = last_hit;
      MOTIF_T *m = h->motif;
      int ws = (translate_dna ? 3 * m->width : m->width);
      //the plus 1 is needed because pos counts from 1
      printed = snprintf(tmp+pos, remaining, "_%ld", length - h->pos - ws + 1);
      if (printed >= remaining) {
        if (pos == 0) die("make_diagram: tmp buffer too small\n");
        break;
      }
      pos += printed;
      remaining -= printed;
      i++;
    }
    for (; i < count; i += (mode == 1 ? 1 : 0), mode *= -1) {
      HIT_T *h = hits[i];
      MOTIF_T *m;
      //skip hits on the wrong strand
      while (h->strand * strand < 0) {
        if (++i >= count) break;
        h = hits[i];
      }
      if (i >= count) break;
      //so we have a hit on this strand
      m = h->motif;
      if (mode == -1) {
        if (h->gap == 0) continue;
        printed = snprintf(tmp+pos, remaining, "_%ld", h->gap);
      } else {
        char *strand = (strand_handling == Combine ? (h->strand == 1 ? "+" : "-") : "");
        int frame = translate_dna ? h->pos % 3 + 1 : 0;
        int offset = 0;
        if (h->pos > 1) {//remember pos is indexed from 1
          if (remaining <= 1) break;
          tmp[pos] = '_';
          offset = 1;
        }
        printed = make_block(m->num, strand, frame, max_hit_pvalue, h->pvalue, print_p, tmp+(pos+offset), remaining - offset) + offset;
      }
      if (printed >= remaining) {//didn't fit
        //error check
        if (pos == 0) die("make_diagram: tmp buffer too small\n");
        break;
      }
      pos += printed;
      remaining -= printed;
      last_hit = h;
    }
    //copy contents of tmp to diagram
    diagram = mm_realloc(diagram, sizeof(char) * (size + pos + 1));
    memcpy(diagram+size, tmp, sizeof(char) * pos);
    size += pos;
    diagram[size] = '\0';
    pos = 0;
    remaining = DIAG_TEMP_SIZE;
  } while (i <= count && strand_hits != 0);

  return diagram;
}

void print_right_trimmed_string(FILE *mast_out, char *text, int len) {
  int i;
  for (i = len-1; i >= 0; --i) {
    if (!isspace(text[i])) break;
  }
  if (i == -1) return;
  fprintf(mast_out,"%.*s",i+1, text); 
}

void append_to_tables(PARSER_STATE_T *ps, SEQUENCE_T *seq, SCORE_T *score) {
  int kmotifs = 0; //disabled

  char *diagram = make_diagram(ps->strand_handling, ps->translate_dna, ps->max_hit_pvalue, FALSE, seq->length, score->strand, seq->hits, seq->hit_pos);

  /* 
    Print the hit in the hits section. 
  */
  if (ps->hit_file != NULL) {
    char *kp = !kmotifs ? "" : (seq->known ? "+ " : "- ");  /* +/- */
    int nl = !kmotifs ? MMSN : MMSN - 2;		/* length of name */
    char *frame = ps->translate_dna ? (score->frame == 0 ? "a" : (score->frame == 1 ? "b" : "c")) : "";
    char *st = ps->strand_handling==Separate ? (score->strand==-1 ? " -" : " +") : ""; 
    int il = (*st=='\0') ? MAXID : MAXID-2;		/* length of id */
    char *elipsis = strlen(seq->comment) > il ? "... " : "   ";
    if (*frame!='\0') il -= 2;
    fprintf(ps->hit_file, "%s%-*.*s %-*.*s%3s%s%s  %8.2g %6ld\n", kp, nl, nl, 
      seq->name, il, il, seq->comment, elipsis, st, frame, score->evalue, seq->length);
  } /* print hit */

  /* 
    Print the motif diagram in the diagrams section.
  */
  if (ps->diag_file != NULL) {
    char hdr[80];						/* header */
    char *st = ps->strand_handling != Separate ? "" : (score->strand == -1 ? " -" : " +"); 
    sprintf(hdr, "%-*.*s%s %8.2g  ", MMSN, MMSN, seq->name, st, score->evalue);
    print_diagram(diagram, hdr, ps->diag_file);
  }

  if (ps->note_file != NULL) {
    char *s;
    int i, gb, commentlines;
    char *sstr = ps->strand_handling != Separate ? "" : (score->strand == -1 ? " (- strand)" : " (+ strand)");
    fprintf(ps->note_file, "\n\n%s%s\n", seq->name, sstr);
    commentlines = 0; 
    for (s = seq->comment, gb=0; ; s+=(gb+1), gb=0) {		/* more to print */
      /* find the next good breaking point */
      for (i = 0; i < PAGEWIDTH - 2; i++) {
        if (s[i] == ' ' || s[i] == '\0') gb = i;		/* good break: blank or EOS */
        if (s[i] == '\0') break; 			/* end of sequence reached */
      }
      if (gb == 0) gb = i - 1;			/* no good break; back up 1 */
      fprintf(ps->note_file, "  %-.*s\n", gb + 1, s);	/* print to break inclusive */ 
      if (s[i] == '\0') break;				/* done with sequence */
      if (++commentlines == MAXIDLINES) {
        fprintf(ps->note_file,"  ...\n"); 
        break;
      }
    }
    
    /* 
      print the length, combined p-value and E-value 
    */
    fprintf(ps->note_file, 
      "  LENGTH = %ld  COMBINED P-VALUE = %8.2e  E-VALUE = %8.2g\n",
      seq->length, score->combined_pvalue, score->evalue);

    /* 
      print the motif diagram 
    */
    print_diagram(diagram, "  DIAGRAM: ", ps->note_file);
    putc('\n', ps->note_file);			/* blank line */

    for (i = 0; i < seq->seg_pos; ++i) {
      int j, k, l, *use_chunk;
      char *num_line, *pvalue_line, *consensus_line, *match_line, *translated_line, *seq_line;
      SEG_T *seg = seq->segs+i;
      int line_len;
      long pos;
      
      seq_line = seg->data;
      line_len = strlen(seq_line);

      num_line = mm_malloc(sizeof(char) * (line_len+1));
      pvalue_line = mm_malloc(sizeof(char) * (line_len+1));
      consensus_line = mm_malloc(sizeof(char) * (line_len+1));
      match_line = mm_malloc(sizeof(char) * (line_len+1));
      if (ps->translate_dna) {
        translated_line = mm_malloc(sizeof(char) * (line_len+1));
      } else {
        translated_line = NULL;
      }
      k = ((int)(line_len / (PAGEWIDTH - 5)) + 1);
      use_chunk = mm_malloc(sizeof(int) * k);
      for (j = 0; j < k; ++j) use_chunk[j] = 0;

      //set all the lines to spaces
      for (j = 0; j < line_len; ++j) {
        num_line[j] = ' ';
        pvalue_line[j] = ' ';
        consensus_line[j] = ' ';
        match_line[j] = ' ';
        if (translated_line) translated_line[j] = ' ';
      }
      num_line[line_len] = '\0';
      pvalue_line[line_len] = '\0';
      consensus_line[line_len] = '\0';
      match_line[line_len] = '\0';
      if (translated_line) translated_line[line_len] = '\0';

      for (j = 0; j < seg->count; ++j) {
        HIT_T *h = seg->hits[j];
        MOTIF_T *m = h->motif;
        int offset = h->pos - seg->start;
        int remaining = line_len - offset + 1;
        int printed;

        //check the hit is on this strand
        if (h->strand * score->strand < 0) continue;

        use_chunk[(int)(offset / (PAGEWIDTH - 5))] += 1;

        //print the motif num
        {
          char *strand = (ps->strand_handling == Combine ? (h->strand == 1 ? "+" : "-") : "");
          int frame = ps->translate_dna ? h->pos % 3 + 1 : 0;
          printed = make_block(m->num, strand, frame, ps->max_hit_pvalue, h->pvalue, FALSE, num_line+offset, remaining);
          //get rid of the poisonous null byte as long as it's not the last byte
          if (printed < (remaining-1)) num_line[offset+printed] = ' ';
        }
        //print the pvalue
        printed = snprintf(pvalue_line+offset, remaining, "%.1e", h->pvalue);
        //get rid of the poisonous null byte as long as it's not the last byte
        if (printed < (remaining-1)) pvalue_line[offset+printed] = ' ';

        //print the consensus, match and translation if applicable
        {
          int end = offset + m->ws;
          int inc = (ps->translate_dna ? 3 : 1);
          char *consensus = (h->strand == 1 ? m->best_f : m->best_r);
          for (k = offset, l = 0; k < end; k +=inc, ++l) {
            consensus_line[k] = consensus[l];
            match_line[k] = h->match[l];
            if (ps->translate_dna) {
              consensus_line[k+1] = '.';
              consensus_line[k+2] = '.';
              translated_line[k] = h->translated[l];
              translated_line[k+1] = '.';
              translated_line[k+2] = '.';
            }
          }
        }
      }
      //output the lines
      int line_width = (PAGEWIDTH - 5);
      for (j = 0, k = 0; j < line_len; j += line_width) {
        line_width = MIN(line_len - j, line_width);
        if (!(use_chunk[(int)(j / (PAGEWIDTH - 5))] || consensus_line[j] != ' ' || consensus_line[j+line_width-1] != ' ')) continue;
        if ((seg->start + j) != 1 || k != 0) putc('\n', ps->note_file);			/* blank line */
        //find the first non-space to copy to the output
        fprintf(ps->note_file, "     ");
        print_right_trimmed_string(ps->note_file, num_line+j, line_width);
        putc('\n', ps->note_file);
        fprintf(ps->note_file, "     ");
        print_right_trimmed_string(ps->note_file, pvalue_line+j, line_width);
        putc('\n', ps->note_file);
        fprintf(ps->note_file, "     ");
        print_right_trimmed_string(ps->note_file, consensus_line+j, line_width);
        putc('\n', ps->note_file);
        fprintf(ps->note_file, "     ");
        print_right_trimmed_string(ps->note_file, match_line+j, line_width);
        putc('\n', ps->note_file);
        if (ps->translate_dna) {
          fprintf(ps->note_file, "     ");
          print_right_trimmed_string(ps->note_file, translated_line+j, line_width);
          putc('\n', ps->note_file);
        }
        fprintf(ps->note_file, "%-5ld", seg->start + j);
        fprintf(ps->note_file, "%.*s\n", line_width, seq_line+j);
        ++k;
      }
      free(num_line);
      free(pvalue_line);
      free(consensus_line);
      free(match_line);
      if (translated_line) free(translated_line);
    }
  }
}
/* }}} */

/* sequence temporary storage {{{ */

void ec_fwrite(const void *ptr, size_t ele, size_t num, FILE *file) {
  size_t written;
  written = fwrite(ptr, ele, num, file);
  if (written != num) {
    if (ferror(file)) die("Error occurred while writing %d elements of size %d, error given as: %s\n", num, ele, strerror(ferror(file)));
    die("Unexplained short write while writing %d elements of size %d\n", num, ele);
  }
}

void ec_fread(void *ptr, size_t ele, size_t num, FILE *file) {
  size_t read;
  read = fread(ptr, ele, num, file);
  if (read != num) {
    if (feof(file)) die("Reached end of file while trying to read %d elements of size %d\n", num, ele);
    if (ferror(file)) die("Error occurred while reading %d elements of size %d, error given as: %s\n", num, ele, strerror(ferror(file)));
    die("Unexplained short read while reading %d elements of size %d\n", num, ele);
  }
}

void write_string(FILE *save, char *str) {
  int len;
  if (str != NULL) { 
    len = strlen(str);
    ec_fwrite(&len, sizeof(int), 1, save);
    if (len > 0) ec_fwrite(str, sizeof(char), len, save);
  } else {
    len = -1;
    ec_fwrite(&len, sizeof(int), 1, save);
  }
}

char* read_string(FILE *save) {
  int len;
  char *str;
  ec_fread(&len, sizeof(int), 1, save);
  if (len >= 0) {
    str = (char*)mm_malloc(sizeof(char) * (len + 1));
    if (len > 0) ec_fread(str, sizeof(char), len, save);
    str[len] = '\0';
    return str;
  } else {
    return NULL;
  }
}

void write_double(FILE *save, double dbl) {
  ec_fwrite(&dbl, sizeof(double), 1, save);
}

double read_double(FILE *save) {
  double dbl;
  ec_fread(&dbl, sizeof(double), 1, save);
  return dbl;
}

void write_int(FILE *save, int integer) {
  ec_fwrite(&integer, sizeof(int), 1, save);
}

int read_int(FILE *save) {
  int integer;
  ec_fread(&integer, sizeof(int), 1, save);
  return integer;
}

void write_long(FILE *save, long lng) {
  ec_fwrite(&lng, sizeof(long), 1, save);
}

long read_long(FILE *save) {
  long lng;
  ec_fread(&lng, sizeof(long), 1, save);
  return lng;
}

void write_hit(FILE *save, HIT_T *hit) {
  write_long(save, hit->pos);
  write_long(save, hit->gap);
  write_int(save, hit->strand);
  write_string(save, hit->motif->id);
  write_double(save, hit->pvalue);
  write_string(save, hit->match);
  write_string(save, hit->translated);
}

HIT_T* read_hit(FILE *save, RBTREE_T *motiflookup) {
  HIT_T *hit;
  char *motifid;
  hit = (HIT_T*)mm_malloc(sizeof(HIT_T));
  hit->pos = read_long(save);
  hit->gap = read_long(save);
  hit->strand = read_int(save);
  motifid = read_string(save);
  hit->motif = rbtree_get(motiflookup, motifid);
  hit->pvalue = read_double(save);
  hit->match = read_string(save);
  hit->translated = read_string(save);
  return hit;
}

void write_seg(FILE *save, SEG_T *seg) {
  int i;
  write_long(save, seg->start);
  write_string(save, seg->data);
  write_int(save, seg->length);
  write_int(save, seg->count);
  for (i = 0; i < seg->count; ++i) {
    write_hit(save, seg->hits[i]);
  }
}

void read_seg(FILE *save, SEG_T *seg, RBTREE_T *motiflookup) {
  int i;
  seg->start = read_long(save);
  seg->data = read_string(save);
  seg->length = read_int(save);
  seg->count = read_int(save);
  seg->hits = (HIT_T**)mm_malloc(sizeof(HIT_T*) * seg->count);
  for (i = 0; i < seg->count; ++i) {
    seg->hits[i] = read_hit(save, motiflookup);
  }
}

long store_sequence(FILE *save, SEQUENCE_T *seq) {
  int i;
  long file_loc;
  file_loc = ftell(save);
  write_string(save, seq->id);
  write_string(save, seq->db->id);
  write_string(save, seq->name);
  write_string(save, seq->comment);
  write_long(save, seq->length);
  write_int(save, seq->known);
  write_int(save, seq->seg_pos);
  for (i = 0; i < seq->seg_pos; ++i) {
    write_seg(save, seq->segs+i);
  }
  return file_loc;
}

SEQUENCE_T* retrieve_sequence(FILE *save, long fp, RBTREE_T *dblookup, RBTREE_T *motiflookup) {
  SEQUENCE_T *seq;
  char *db_id;
  int i, j, k, hit_count;
  long where; //where we are now
  seq = mm_malloc(sizeof(SEQUENCE_T));
  memset(seq, 0, sizeof(SEQUENCE_T));
  where = ftell(save);
  fseek(save, fp, SEEK_SET);
  seq->id = read_string(save);
  db_id = read_string(save);
  seq->db = rbtree_get(dblookup, db_id);
  free(db_id);
  seq->name = read_string(save);
  seq->comment = read_string(save);
  seq->length = read_long(save);
  seq->known = read_int(save);
  seq->seg_pos = read_int(save);
  seq->seg_count = seq->seg_pos;
  seq->segs = (SEG_T*)mm_malloc(sizeof(SEG_T) * seq->seg_count);
  for (i = 0, hit_count = 0; i < seq->seg_pos; ++i) {
    SEG_T *seg = seq->segs+i;
    read_seg(save, seg, motiflookup);
    hit_count += seg->count;
  }
  seq->hit_count = hit_count;
  seq->hit_pos = hit_count;
  seq->hits = (HIT_T**)mm_malloc(sizeof(HIT_T*) * hit_count);
  k = 0;
  for (i = 0; i < seq->seg_pos; ++i) {
    SEG_T *seg = seq->segs+i;
    for (j = 0; j < seg->count; ++j) {
      seq->hits[k++] = seg->hits[j];
    }
  }
  //before returning move the file pointer back to its original location
  fseek(save, where, SEEK_SET);
  return seq;
}

void destroy_sequence(SEQUENCE_T *seq) {
  int i;
  free(seq->id);
  free(seq->name);
  free(seq->comment);
  for (i = 0; i < seq->seg_pos; ++i) {
    SEG_T *seg = seq->segs+i;
    free(seg->data);
    free(seg->hits);
  }
  free(seq->segs);
  for (i = 0; i < seq->hit_pos; ++i) {
    HIT_T *hit = seq->hits[i];
    free(hit->match);
    free(hit->translated);
    free(hit);
  }
  free(seq->hits);
}

/* }}} */

/* parsing helper functions {{{*/
int compare_strings(void *v1, void *v2) {
  char *str1, *str2;
  str1 = (char*)v1;
  str2 = (char*)v2;
  return strcmp(str1, str2);
}

int compare_pstrings(const void *v1, const void *v2) {
  char *str1, *str2;
  str1 = *((char**)v1);
  str2 = *((char**)v2);
  return strcmp(str1, str2);
}

int compare_scores(void *v1, void *v2) {
  SCORE_T *sc1, *sc2;
  double diff;
  int db_diff, seq_diff;
  long fp_diff;

  sc1 = (SCORE_T*)v1;
  sc2 = (SCORE_T*)v2;
  //order by score
  diff = sc1->combined_pvalue - sc2->combined_pvalue;
  if (diff < 0) {
    return -1;
  } else if (diff > 0) {
    return 1;
  }
  //order by database
  db_diff = sc1->db_num - sc2->db_num;
  if (db_diff < 0) {
    return -1;
  } else if (db_diff > 0) {
    return 1;
  }
  //order by db sequence number
  seq_diff = sc1->seq_num - sc2->seq_num;
  if (seq_diff < 0) {
    return -1;
  } else if (seq_diff > 0) {
    return 1;
  }
  //must be the same sequence
  //order by strand
  if (sc1->strand > sc2->strand) {
    return -1;
  } else if (sc1->strand < sc2->strand) {
    return 1;
  }
  DEBUG_MSG(QUIET_VERBOSE, "Warning, compare_scores compared two scores as identical. This is probably an error.\n");
  //I don't think this is possible
  return 0;
}

int ld_str(char *value, void *data) {
  copy_string((char**)data, value);
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
  parsed_value = strtol(value, &end_ptr, 10);
  if (end_ptr == value) return -1;
  if (errno) return errno;
  *((long*)data) = parsed_value;
  return 0;
}

int ld_double(char *value, void *data) {
  double parsed_value;
  char *end_ptr;
  parsed_value = strtod(value, &end_ptr);
  if (end_ptr == value) return -1;
  if (errno) return errno;
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

void parse_attributes(PARSER_STATE_T *ps, char *tag, const xmlChar **attrs, 
    int count, char **names, int (**parsers)(char*, void*), void **parser_data, BOOLEAN_T *required) {
  int i, found;
  BOOLEAN_T done[count];
  for (i = 0; i < count; ++i) done[i] = FALSE;
  for (i = 0; attrs[i] != NULL; i += 2) {
    found = binary_search(attrs+i, names, count, sizeof(char*), compare_pstrings);
    if (found >= 0) {
      int (*parser)(char*, void*);
      void* data;
      if (done[found]) {
        ps->state = PS_ERROR;
        DEBUG_FMT(QUIET_VERBOSE, "Already parsed %s::%s.\n", tag, names[found]);
        continue;
      }
      done[found] = TRUE;
      parser = parsers[found];
      data = parser_data[found];
      if (parser((char*)attrs[i+1], data)) {
        ps->state = PS_ERROR;
        DEBUG_FMT(QUIET_VERBOSE, "Bad value \"%s\" for %s::%s.\n", (char*)attrs[i+1], tag, names[found]);
        continue;
      }
    }
  }
  for (i = 0; i < count; ++i) {
    if (required[i] && !done[i]) {
      ps->state = PS_ERROR;
      DEBUG_FMT(QUIET_VERBOSE, "Missing required attribute %s::%s.\n", tag, names[i]);
    }
  }
}
/*}}}*/

/* parse mast start tag {{{*/
void start_ele_mast(PARSER_STATE_T *ps, const xmlChar **attrs) {
  char* names[2] = {"release", "version"};
  int (*parsers[2])(char*, void*) = {ld_str, ld_str};
  void *data[2] = {&(ps->release), &(ps->version)};
  BOOLEAN_T required[2] = {TRUE, TRUE};
  parse_attributes(ps, "mast", attrs, 2, names, parsers, data, required);
}
/*}}}*/

/* parse model {{{*/
void end_ele_command_line(PARSER_STATE_T *ps) {
  copy_string(&(ps->command_line), ps->characters.buffer);
}

void end_ele_max_correlation(PARSER_STATE_T *ps) {
  if (ld_double(ps->characters.buffer, &(ps->max_correlation))) {
    ps->state = PS_ERROR;
    DEBUG_FMT(QUIET_VERBOSE, "Error parsing max correlation value \"%s\" as double.\n", ps->characters.buffer);
  }
}

void start_ele_remove_correlated(PARSER_STATE_T *ps, const xmlChar **attrs) {
  char* opts_value[2] = {"n", "y"};
  int outs_value[2] = {FALSE, TRUE};
  MULTI_T multi_value = {.count = 2, .options = opts_value, .outputs = outs_value, .target = &(ps->remove_correlated)};
  char* names[1] = {"value"};
  int (*parsers[1])(char*, void*) = {ld_multi};
  void *data[1] = {&multi_value};
  BOOLEAN_T required[1] = {TRUE};
  parse_attributes(ps, "remove_correlated", attrs, 1, names, parsers, data, required);
}

void start_ele_strand_handling(PARSER_STATE_T *ps, const xmlChar **attrs) {
  char* opts_value[4] = {"combine", "norc", "protein", "separate"};
  int outs_value[4] = {Combine, Norc, Protein, Separate};
  MULTI_T multi_value = {.count = 4, .options = opts_value, .outputs = outs_value, .target = &(ps->strand_handling)};
  char* names[1] = {"value"};
  int (*parsers[1])(char*, void*) = {ld_multi};
  void *data[1] = {&multi_value};
  BOOLEAN_T required[1] = {TRUE};
  parse_attributes(ps, "strand_handling", attrs, 1, names, parsers, data, required);
}

void start_ele_translate_dna(PARSER_STATE_T *ps, const xmlChar **attrs) {
  char* opts_value[2] = {"n", "y"};
  int outs_value[2] = {FALSE, TRUE};
  MULTI_T multi_value = {.count = 2, .options = opts_value, .outputs = outs_value, .target = &(ps->translate_dna)};
  char* names[1] = {"value"};
  int (*parsers[1])(char*, void*) = {ld_multi};
  void *data[1] = {&multi_value};
  BOOLEAN_T required[1] = {TRUE};
  parse_attributes(ps, "translate_dna", attrs, 1, names, parsers, data, required);
}

void end_ele_max_seq_evalue(PARSER_STATE_T *ps) {
  if (ld_double(ps->characters.buffer, &(ps->max_seq_evalue))) {
    ps->state = PS_ERROR;
    DEBUG_FMT(QUIET_VERBOSE, "Error parsing max sequence evalue \"%s\" as double.\n", ps->characters.buffer);
  }
}

void start_ele_adj_hit_pvalue(PARSER_STATE_T *ps, const xmlChar **attrs) {
  char* opts_value[2] = {"n", "y"};
  int outs_value[2] = {FALSE, TRUE};
  MULTI_T multi_value = {.count = 2, .options = opts_value, .outputs = outs_value, .target = &(ps->adj_hit_pvalue)};
  char* names[1] = {"value"};
  int (*parsers[1])(char*, void*) = {ld_multi};
  void *data[1] = {&multi_value};
  BOOLEAN_T required[1] = {TRUE};
  parse_attributes(ps, "adj_hit_pvalue", attrs, 1, names, parsers, data, required);
}

void end_ele_max_hit_pvalue(PARSER_STATE_T *ps) {
  if (ld_double(ps->characters.buffer, &(ps->max_hit_pvalue))) {
    ps->state = PS_ERROR;
    DEBUG_FMT(QUIET_VERBOSE, "Error parsing max hit pvalue \"%s\" as double.\n", ps->characters.buffer);
  }
}

void end_ele_max_weak_pvalue(PARSER_STATE_T *ps) {
  if (ld_double(ps->characters.buffer, &(ps->max_weak_pvalue))) {
    ps->state = PS_ERROR;
    DEBUG_FMT(QUIET_VERBOSE, "Error parsing max weak-hit pvalue \"%s\" as double.\n", ps->characters.buffer);
  }
}

void end_ele_host(PARSER_STATE_T *ps) {
  copy_string(&(ps->host), ps->characters.buffer);
}

void end_ele_when(PARSER_STATE_T *ps) {
  copy_string(&(ps->when), ps->characters.buffer);
}

/*}}}*/

/* parse alphabet {{{*/
void start_ele_alphabet(PARSER_STATE_T *ps, const xmlChar **attrs) {
  char* opts_bg_source[3] = {"file", "preset", "sequence_composition"};
  int outs_bg_source[3] = {BGFILE, PRESET, SEQCOMP};
  MULTI_T multi_bg_source = {.count = 3, .options = opts_bg_source, .outputs = outs_bg_source, .target = &(ps->alphabet.source)};
  char* opts_type[2] = {"amino-acid", "nucleotide"};
  int outs_type[2] = {FALSE, TRUE};
  MULTI_T multi_type = {.count = 2, .options = opts_type, .outputs = outs_type, .target = &(ps->alphabet.dna)};
  char* names[3] = {"bg_file", "bg_source", "type"};
  int (*parsers[3])(char*, void*) = {ld_str, ld_multi, ld_multi};
  void *data[3] = {&(ps->alphabet.bg_file), &multi_bg_source, &multi_type};
  BOOLEAN_T required[3] = {FALSE, TRUE, TRUE};
  parse_attributes(ps, "alphabet", attrs, 3, names, parsers, data, required);
  if (ps->state == PS_IN_ALPHABET) {
    if (ps->translate_dna) {
      if (!(ps->alphabet.dna)) {
        ps->db_dna = TRUE;
      } else {
        ps->state = PS_ERROR;
        DEBUG_MSG(QUIET_VERBOSE, "Alphabet type must be protein when translating dna.\n");
      }
    } else {
      ps->db_dna = ps->alphabet.dna;
    }
    if (ps->alphabet.source == BGFILE && ps->alphabet.bg_file == NULL) {
      ps->state = PS_ERROR;
      DEBUG_MSG(QUIET_VERBOSE, "Alphabet bg_source is file but no bg_file attribute is specified.\n");
    }
  }
}

void start_ele_letter(PARSER_STATE_T *ps, const xmlChar **attrs) {
  ALPHABET_T *alph = &(ps->alphabet);
  ALPHSYM_T *sym;
  alph->count += 1;
  alph->symbols = mm_realloc(alph->symbols, sizeof(ALPHSYM_T) * alph->count);
  sym = (alph->symbols)+(alph->count - 1);

  char* opts_ambig[2] = {"n", "y"};
  int outs_ambig[2] = {FALSE, TRUE};
  MULTI_T multi_ambig = {.count = 2, .options = opts_ambig, .outputs = outs_ambig, .target = &(sym->ambig)};

  char* names[3] = {"ambig", "bg_value", "symbol"};
  int (*parsers[3])(char*, void*) = {ld_multi, ld_double, ld_char};
  void *data[3] = {&multi_ambig, &(sym->bg_value), &(sym->symbol)};
  BOOLEAN_T required[3] = {FALSE, FALSE, TRUE};
  parse_attributes(ps, "letter", attrs, 3, names, parsers, data, required);
}
/*}}}*/

/* parse motifs {{{*/
void start_ele_motifs(PARSER_STATE_T *ps, const xmlChar **attrs) {
  char* names[3] = {"last_mod_date", "name", "source"};
  int (*parsers[3])(char*, void*) = {ld_str, ld_str, ld_str};
  void *data[3] = {&(ps->motifs.last_mod_date), &(ps->motifs.name), &(ps->motifs.source)};
  BOOLEAN_T required[3] = {TRUE, TRUE, TRUE};
  parse_attributes(ps, "motifs", attrs, 3, names, parsers, data, required);
}

void start_ele_motif(PARSER_STATE_T *ps, const xmlChar **attrs) {
  MOTIFS_T *motifs = &(ps->motifs);
  MOTIF_T *motif;
  int len;
  
  motif = (MOTIF_T*)mm_malloc(sizeof(MOTIF_T));
  memset(motif, 0, sizeof(MOTIF_T)); 
  motifs->count += 1;
  motifs->all_motifs = mm_realloc(motifs->all_motifs, sizeof(MOTIF_T*) * motifs->count);
  motifs->all_motifs[motifs->count - 1] = motif;

  char* opts_bad[2] = {"n", "y"};
  int outs_bad[2] = {FALSE, TRUE};
  MULTI_T multi_bad = {.count = 2, .options = opts_bad, .outputs = outs_bad, .target = &(motif->bad)};

  char* names[7] = {"bad", "best_f", "best_r", "id", "name", "num", "width"}; 
  int (*parsers[7])(char*, void*) = {ld_multi, ld_str, 
    ld_str, ld_str, ld_str, 
    ld_int, ld_int};
  void *data[7] = {&multi_bad, &(motif->best_f), &(motif->best_r), &(motif->id), &(motif->name), &(motif->num), &(motif->width)};
  BOOLEAN_T required[7] = {FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE};
  parse_attributes(ps, "motif", attrs, 7, names, parsers, data, required);

  if (ps->state == PS_IN_MOTIF) {
    //check the length of the consensus
    len = strlen(motif->best_f);
    if (len != motif->width) {
      ps->state = PS_ERROR;
      DEBUG_FMT(QUIET_VERBOSE, "motif::best_f \"%s\" is a different length to the motif (%d)\n", motif->best_f, motif->width);
      return;
    }
    if (motif->best_r) {
      len = strlen(motif->best_r);
      if (len != motif->width) {
        ps->state = PS_ERROR;
        DEBUG_FMT(QUIET_VERBOSE, "motif::best_r \"%s\" is a different length to the motif (%d)\n", motif->best_r, motif->width);
        return;
      }
    }

    if (ps->translate_dna) {
      int i, j;
      motif->ws = motif->width * 3;

      //use the reversed best_f string for the best_r
      motif->best_r = mm_malloc(sizeof(char) * (motif->width + 1));
      for (i = 0, j = motif->width-1; i < motif->width; ++i, --j) {
        motif->best_r[i] = motif->best_f[j];
      }
      motif->best_r[motif->width] = '\0';
    } else {
      motif->ws = motif->width;
    }

    if (motif->num > 1) {
      motif->correlations = mm_malloc(sizeof(double) * (motif->num - 1));
    }

    rbtree_put(motifs->lookup, motif->id, motif);
  }
}

void start_ele_correlation(PARSER_STATE_T *ps, const xmlChar **attrs) {
  char *id_a, *id_b;
  double correlation;
  MOTIF_T *motif_a, *motif_b;
  id_a = NULL;
  id_b = NULL;
  char* names[3] = {"motif_a", "motif_b", "value"};
  int (*parsers[3])(char*, void*) = {ld_str, ld_str, ld_double};
  void *data[3] = {&id_a, &id_b, &correlation};
  BOOLEAN_T required[3] = {TRUE, TRUE, TRUE};
  parse_attributes(ps, "correlation", attrs, 3, names, parsers, data, required);
  //error check
  if (ps->state == PS_IN_CORRELATION) {
    //find the motifs
    motif_a = rbtree_get(ps->motifs.lookup, id_a);
    motif_b = rbtree_get(ps->motifs.lookup, id_b);
    if (motif_a && motif_b) {
      if (motif_a->num > motif_b->num) {
        motif_a->correlations[motif_b->num - 1] = correlation;
      } else {
        motif_b->correlations[motif_a->num - 1] = correlation;
      }
    } else {
      ps->state = PS_ERROR;
      if (motif_a == NULL) {
        DEBUG_FMT(QUIET_VERBOSE, "correlation::motif_a \"%s\" is a bad reference\n", id_a);
      }
      if (motif_b == NULL) {
        DEBUG_FMT(QUIET_VERBOSE, "correlation::motif_b \"%s\" is a bad reference\n", id_b); 
      }
    }
  }
  //clean up
  if (id_a) free(id_a);
  if (id_b) free(id_b);
}

void start_ele_nos(PARSER_STATE_T *ps, const xmlChar **attrs) {
  char *names[1] = {"length"};
  int (*parsers[1])(char*, void*) = {ld_int};
  void *data[1] = {&(ps->motifs.diagram_remain)};
  BOOLEAN_T required[1] = {TRUE};
  parse_attributes(ps, "nos", attrs, 1, names, parsers, data, required);
}

void start_ele_expect(PARSER_STATE_T *ps, const xmlChar **attrs) {
  int pos, gap;
  char *motif_id = NULL;
  MOTIF_T *motif;
  char *names[3] = {"gap", "motif", "pos"};
  int (*parsers[3])(char*, void*) = {ld_int, ld_str, ld_int};
  void *data[3] = {&gap, &motif_id, &pos};
  BOOLEAN_T required[3] = {TRUE, TRUE, TRUE};
  parse_attributes(ps, "expect", attrs, 3, names, parsers, data, required);
  if (ps->state == PS_IN_EXPECT) {
    motif = rbtree_get(ps->motifs.lookup, motif_id);
    if (motif) {
      char *diagram;
      char tmp[40];
      int printed;
      if (ps->motifs.diagram_motifs == 0) {
        if (gap == 0) {
          printed = snprintf(tmp, 40, "[%d]", motif->num);
        } else {
          printed = snprintf(tmp, 40, "%d_[%d]", gap, motif->num);
        }
      } else {
        if (gap == 0) {
          printed = snprintf(tmp, 40, "_[%d]", motif->num);
        } else {
          printed = snprintf(tmp, 40, "_%d_[%d]", gap, motif->num);
        }
      }
      diagram = ps->motifs.diagram;
      if (diagram) {
        int len;
        len = strlen(diagram);
        diagram = mm_realloc(diagram, sizeof(char) * (len + printed + 1));
        strcpy(diagram+len, tmp);
      } else {
        diagram = mm_malloc(sizeof(char) * (printed + 1));
        strcpy(diagram, tmp);
      }
      ps->motifs.diagram_motifs += 1;
      ps->motifs.diagram_remain -= gap;
      ps->motifs.diagram_remain -= motif->ws;
    } else {
      ps->state = PS_ERROR;
      DEBUG_FMT(QUIET_VERBOSE, "expect::motif \"%s\" is a bad reference\n", motif_id);
    }
  }
  if (motif_id != NULL) free(motif_id);
}

void end_ele_nos(PARSER_STATE_T *ps) {
  char *diagram = ps->motifs.diagram;
  int gap = ps->motifs.diagram_remain;
  if (gap > 0) {
    char tmp[40];
    int printed;
    printed = snprintf(tmp, 40, "_%d", gap);
    if (diagram) {
      int len;
      len = strlen(diagram);
      diagram = mm_realloc(diagram, sizeof(char) * (len + printed + 1));
      strcpy(diagram+len, tmp);
    } else {
      diagram = mm_malloc(sizeof(char) * (printed + 1));
      strcpy(diagram, tmp);
    }
  } else if (gap < 0) {
    ps->state = PS_ERROR;
    DEBUG_MSG(QUIET_VERBOSE, "The nominal order and spacing diagram extends past the length.\n");
  }
}

/*}}}*/

/* parse sequences {{{*/
void start_ele_database(PARSER_STATE_T *ps, const xmlChar **attrs) {
  DATABASES_T *databases = &(ps->databases);
  DATABASE_T *db;
  db = mm_malloc(sizeof(DATABASE_T));
  memset(db, 0, sizeof(DATABASE_T));
  databases->count += 1;
  databases->all_databases = mm_realloc(databases->all_databases, sizeof(DATABASE_T*) * databases->count);
  databases->all_databases[databases->count - 1] = db;

  char* opts_type[2] = {"amino-acid", "nucleotide"};
  int outs_type[2] = {FALSE, TRUE};
  MULTI_T multi_type = {.count = 2, .options = opts_type, .outputs = outs_type, .target = &(db->is_dna)};

  char *names[8] = {"id", "last_mod_date", "name", "num", "residue_count", "seq_count", "source", "type"};
  int (*parsers[8])(char*, void*) = {ld_str, ld_str, ld_str, ld_int, ld_int, ld_int, ld_str, ld_multi};
  void *data[8] = {&(db->id), &(db->last_mod_date), &(db->name), &(db->num), &(db->residue_count), &(db->seq_count), &(db->source), &multi_type};
  BOOLEAN_T required[8] = {TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE};
  parse_attributes(ps, "database", attrs, 8, names, parsers, data, required);
  rbtree_put(databases->lookup, db->id, db);
}

void start_ele_sequence(PARSER_STATE_T *ps, const xmlChar **attrs) {
  SEQUENCE_T *seq = &(ps->current_seq);
  char *db_id = NULL;
  
  //seq defaults, we can't just zero everything with memset because of the segs and hits we cache
  seq->id = NULL;
  seq->db = NULL;
  seq->name = NULL;
  seq->num = -1;
  seq->comment = NULL;
  seq->length = 0;
  seq->known = FALSE;
  seq->has_score1 = FALSE;
  seq->has_score2 = FALSE;

  char *names[6] = {"comment","db","id","length","name", "num"};
  int (*parsers[6])(char*, void*) = {ld_str, ld_str, ld_str, ld_long, ld_str, ld_int};
  void *data[6] = {&(seq->comment), &db_id, &(seq->id), &(seq->length), &(seq->name), &(seq->num)};
  BOOLEAN_T required[6] = {FALSE, TRUE, TRUE, TRUE, TRUE, TRUE};
  parse_attributes(ps, "sequence", attrs, 6, names, parsers, data, required);
  if (ps->state == PS_IN_SEQUENCE) {
    seq->db = rbtree_get(ps->databases.lookup, db_id);
    if (!seq->db) {
      ps->state = PS_ERROR;
      DEBUG_FMT(QUIET_VERBOSE, "sequence::db \"%s\" is a bad reference\n", db_id);
    }
  }
  if (db_id) free(db_id);
}

void start_ele_score(PARSER_STATE_T *ps, const xmlChar **attrs) {
  int strand, frame = 0;
  double combined_pvalue, evalue;

  char* opts_frame[3] = {"a", "b", "c"};
  int outs_frame[3] = {0, 1, 2};
  MULTI_T multi_frame = {.count = 3, .options = opts_frame, .outputs = outs_frame, .target = &frame};

  char* opts_strand[3] = {"both", "forward", "reverse"};
  int outs_strand[3] = {0, 1, -1};
  MULTI_T multi_strand = {.count = 3, .options = opts_strand, .outputs = outs_strand, .target = &strand};

  char *names[4] = {"combined_pvalue", "evalue", "frame", "strand"};
  int (*parsers[4])(char*, void*) = {ld_double, ld_double, ld_multi, ld_multi};
  void *data[4] = {&combined_pvalue, &evalue, &multi_frame, &multi_strand};
  BOOLEAN_T required[4] = {TRUE, TRUE, FALSE, TRUE};
  parse_attributes(ps, "score", attrs, 4, names, parsers, data, required);
  if (ps->state == PS_IN_SCORE) {
    SCORE_T *score;
    if (strand == 0 || strand == 1) {
      score = &(ps->current_seq.score1);
      ps->current_seq.has_score1 = TRUE;
    } else {
      score = &(ps->current_seq.score2);
      ps->current_seq.has_score2 = TRUE;
    }
    score->strand = strand;
    score->frame = frame;
    score->combined_pvalue = combined_pvalue;
    score->evalue = evalue;
    score->db_num = ps->current_seq.db->num;
    score->seq_num = ps->current_seq.num;
    score->file_loc = -1;
  }
  ps->databases.total_hit_count += 1;
}

void start_ele_seg(PARSER_STATE_T *ps, const xmlChar **attrs) {
  SEQUENCE_T *seq = &(ps->current_seq);
  SEG_T *seg;
  seq->seg_pos += 1;
  if (seq->seg_count < seq->seg_pos) {
    seq->seg_count = seq->seg_pos;
    seq->segs = mm_realloc(seq->segs, sizeof(SEG_T) * seq->seg_count);
  }
  seg = seq->segs+(seq->seg_pos - 1);
  memset(seg, 0, sizeof(SEG_T));
  
  char *names[1] = {"start"};
  int (*parsers[1])(char*, void*) = {ld_long};
  void *data[1] = {&(seg->start)};
  BOOLEAN_T required[1] = {TRUE};
  parse_attributes(ps, "seg", attrs, 1, names, parsers, data, required);
}

void end_ele_data(PARSER_STATE_T *ps) {
  SEG_T *seg = ps->current_seq.segs+(ps->current_seq.seg_pos - 1);
  copy_string(&(seg->data), ps->characters.buffer);
  seg->length = strlen(seg->data);
}

void start_ele_hit(PARSER_STATE_T *ps, const xmlChar **attrs) {
  SEQUENCE_T *seq = &(ps->current_seq);
  SEG_T *seg = seq->segs+(seq->seg_pos - 1);
  HIT_T *hit;
  seq->hit_pos += 1;
  if (seq->hit_count < seq->hit_pos) {
    hit = (HIT_T*)mm_malloc(sizeof(HIT_T));
    seq->hits = (HIT_T**)mm_realloc(seq->hits, sizeof(HIT_T*) * seq->hit_pos);
    seq->hits[seq->hit_count++] = hit;
  } else {
    hit = seq->hits[seq->hit_pos - 1];
  }
  memset(hit, 0, sizeof(HIT_T));
  
  hit->strand = 1;//default
  char* opts_strand[3] = {"forward", "reverse"};
  int outs_strand[3] = {1, -1};
  MULTI_T multi_strand = {.count = 3, .options = opts_strand, .outputs = outs_strand, .target = &(hit->strand)};

  char *motif_id = NULL;
  char *names[7] = {"gap", "match", "motif", "pos", "pvalue", "strand", "translation"};
  int (*parsers[7])(char*, void*) = {ld_long, ld_str, ld_str, ld_long, ld_double, ld_multi, ld_str};
  void *data[7] = {&(hit->gap), &(hit->match), &(motif_id), &(hit->pos), &(hit->pvalue), &multi_strand, &(hit->translated)};
  BOOLEAN_T required[7] = {TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, ps->translate_dna};
  parse_attributes(ps, "hit", attrs, 7, names, parsers, data, required);
  //check for problems
  if (ps->state == PS_IN_HIT) {
    hit->motif = rbtree_get(ps->motifs.lookup, motif_id);
    if (!(hit->motif)) {
      ps->state = PS_ERROR;
      DEBUG_FMT(QUIET_VERBOSE, "hit::motif \"%s\" is a bad reference\n", motif_id);
    } else if (hit->pos < seg->start) {
      //check that the position is within the bounds of the segment
      ps->state = PS_ERROR;
      DEBUG_FMT(QUIET_VERBOSE, "hit::pos (%ld) is before the start (%ld) of the segment\n", hit->pos, seg->start);
    } else if ((hit->pos - seg->start) + hit->motif->ws > seg->length) {
      //check that the hit ends before the bounds of the segment
      ps->state = PS_ERROR;
      DEBUG_MSG(QUIET_VERBOSE, "hit ends after the end of the segment\n");
    }
    seg->count += 1;
  }
  if (motif_id) free(motif_id);
}

void end_ele_seg(PARSER_STATE_T *ps) {
  int i, j;
  SEQUENCE_T *seq = &(ps->current_seq);
  SEG_T *seg = seq->segs+(seq->seg_pos - 1);
  seg->hits = (HIT_T**)mm_malloc(sizeof(HIT_T*) * seg->count);
  for (i = 0, j = seq->hit_pos - seg->count; i < seg->count; ++i, j++) {
    seg->hits[i] = seq->hits[j];
  }
}

void end_ele_sequence(PARSER_STATE_T *ps) {
  int i;
  SEQUENCE_T *cur_seq, *postponed_seq;
  SCORE_T *score1, *score2, *postponed_score, *good_score, *worse_score;
  RBNODE_T *node;
  cur_seq = &(ps->current_seq);
  score1 = cur_seq->has_score1 ? &(cur_seq->score1) : NULL;
  score2 = cur_seq->has_score2 ? &(cur_seq->score2) : NULL;

  good_score = NULL;
  worse_score = NULL;
  if (score1 && score2) {
    if (compare_scores(score1, score2) < 0) {
      good_score = score1;
      worse_score = score2;
    } else {
      good_score = score2;
      worse_score = score1;
    }
  } else if (score1) {
    good_score = score1;
  } else if (score2) {
    good_score = score2;
  } else {
    ps->state = PS_ERROR;
    DEBUG_FMT(QUIET_VERBOSE, "No scores for sequence \"%s\".\n", cur_seq->name);
  }
  if (ps->state != PS_ERROR) {
    //output any postponed scores that are better
    node = rbtree_first(ps->postponed_scores);
    postponed_score = (node ? rbtree_key(node) : NULL);
    while (postponed_score && compare_scores(postponed_score, good_score) < 0) {
      rbtree_delete(ps->postponed_scores, node, NULL, NULL);
      //retrieve sequence from file and output score
      postponed_seq = retrieve_sequence(ps->save, postponed_score->file_loc, ps->databases.lookup, ps->motifs.lookup);
      append_to_tables(ps, postponed_seq, postponed_score);
      //clean up
      free(postponed_score);
      destroy_sequence(postponed_seq);
      //get next postponed score
      node = rbtree_first(ps->postponed_scores);
      postponed_score = (node ? rbtree_key(node) : NULL);
    }
    //output the best existing score
    append_to_tables(ps, cur_seq, good_score);
    
    if (worse_score) {
      SCORE_T *copy;
      long file_loc;
      //store sequence information into a temporary file so it can be retrieved later
      file_loc = store_sequence(ps->save, cur_seq);
      //make a copy of the score and postpone it
      copy = (SCORE_T*)mm_malloc(sizeof(SCORE_T));
      copy->combined_pvalue = worse_score->combined_pvalue;
      copy->evalue = worse_score->evalue;
      copy->db_num = worse_score->db_num;
      copy->strand = worse_score->strand;
      copy->frame = worse_score->frame;
      copy->file_loc = file_loc;
      rbtree_put(ps->postponed_scores, copy, NULL);
    }
  }
  //clean up the sequence allocations
  for (i = 0; i < cur_seq->seg_pos; i++) {
    SEG_T *seg = cur_seq->segs+i;
    //clean up seg
    if (seg->data) free(seg->data);
    if (seg->hits) free(seg->hits);
  }
  cur_seq->seg_pos = 0;
  for (i = 0; i < cur_seq->hit_pos; i++) {
    HIT_T *hit = cur_seq->hits[i];
    //clean up hit
    if (hit->match) free(hit->match);
    if (hit->translated) free(hit->translated);
  }
  cur_seq->hit_pos = 0;
  if (cur_seq->name) free(cur_seq->name);
  cur_seq->name = NULL;
  if (cur_seq->id) free(cur_seq->id);
  cur_seq->id = NULL;
  if (cur_seq->comment) free(cur_seq->comment);
  cur_seq->comment = NULL;
}

void end_ele_sequences(PARSER_STATE_T *ps) {
  //output any left over sequences
  RBNODE_T *node;
  SEQUENCE_T *postponed_seq;
  SCORE_T *postponed_score;
  for (node = rbtree_first(ps->postponed_scores); node != NULL; node = rbtree_next(node)) {
    postponed_score = rbtree_key(node);
    //retrieve sequence from file and output score
    postponed_seq = retrieve_sequence(ps->save, postponed_score->file_loc, ps->databases.lookup, ps->motifs.lookup);
    append_to_tables(ps, postponed_seq, postponed_score);
    //clean up
    free(postponed_score);
    destroy_sequence(postponed_seq);
  }
  rbtree_destroy(ps->postponed_scores);
  ps->postponed_scores = NULL;
}
/*}}}*/

/* parse runtime {{{*/
void start_ele_runtime(PARSER_STATE_T *ps, const xmlChar **attrs) {
  char* names[1] = {"seconds"};
  int (*parsers[1])(char*, void*) = {ld_double};
  void *data[1] = {&(ps->runtime)};
  BOOLEAN_T required[1] = {TRUE};
  parse_attributes(ps, "runtime", attrs, 1, names, parsers, data, required);
}
/*}}}*/

/* parse mast close tag - begin output {{{*/
void end_ele_mast(PARSER_STATE_T *ps) {
  output_mast_txt(ps);
}
/*}}}*/

/* main, SAX callbacks and state machine {{{*/
void handle_start_doc(void *ctx) {
  char* txt_file_name;
  PARSER_STATE_T *ps = (PARSER_STATE_T*)ctx;
  txt_file_name = ps->txt_file_name;
  memset(ps, 0, sizeof(PARSER_STATE_T));
  ps->txt_file_name = txt_file_name;
  ps->state = PS_START;
  ps->udepth = 0;
  ps->save = tmpfile();
  ps->hit_file = tmpfile();
  ps->diag_file = tmpfile();
  ps->note_file = tmpfile();
  //set up buffer
  ps->characters.buffer = mm_malloc(sizeof(char)*10);
  ps->characters.size = 10;
  ps->characters.pos = 0;
  ps->characters.buffer[0] = '\0';
  //set up alphabet
  ps->alphabet.symbols = NULL;
  //set up motifs
  ps->motifs.lookup = rbtree_create(compare_strings, NULL, NULL, NULL, NULL);
  //set up dbs
  ps->databases.lookup = rbtree_create(compare_strings, NULL, NULL, NULL, NULL);
  //set up postponed scores
  ps->postponed_scores = rbtree_create(compare_scores, NULL, NULL, NULL, NULL);
}

void handle_end_doc(void *ctx) {
  PARSER_STATE_T *ps = (PARSER_STATE_T*)ctx;
  //TODO cleanup
  if (ps->version) free(ps->version);
  if (ps->release) free(ps->release);
  fclose(ps->save);
  fclose(ps->hit_file);
  fclose(ps->diag_file);
  fclose(ps->note_file);
}

int accept_all(char ch) {
  return 1;
}

int accept_not_isspace(char ch) {
  return !isspace(ch);
}

void handle_characters(void *ctx, const xmlChar *ch, int len) {
  CHARBUF_T *buf;
  int i, start, end, accepted;
  int (*accept) (char);
  PARSER_STATE_T *ps = (PARSER_STATE_T*)ctx;
  accept = ps->characters.accept;
  buf = &(ps->characters);
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

#define DO_START_ELE(_expected_,_transition_,_char_accept_) \
  if (strcmp((char*)name, #_expected_ ) == 0) { \
    ps->state = _transition_; \
    start_ele_ ## _expected_ (ps, attrs); \
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
  PARSER_STATE_T *ps = (PARSER_STATE_T*)ctx;
  int known;
  if (ps->state == PS_ERROR) return;
  if (ps->udepth) {//we don't know where we are!
    ps->udepth += 1;
  } else {
    //reset the character buffer to the begining
    ps->characters.pos = 0;
    ps->characters.buffer[0] = '\0';
    known = 1; //assume we can find it
    switch (ps->state) {
      case PS_START:
        DO_START_ELE(mast,PS_IN_MAST, IGNORE);
        known = 0;
        break;
      case PS_IN_MAST:
        CHECK_START_ELE(model, PS_IN_MODEL, IGNORE);
        DO_START_ELE(alphabet, PS_IN_ALPHABET, IGNORE);
        DO_START_ELE(motifs, PS_IN_MOTIFS, IGNORE);
        CHECK_START_ELE(sequences, PS_IN_SEQUENCES, IGNORE);
        DO_START_ELE(runtime, PS_IN_RUNTIME, IGNORE);
        known = 0;
        break;
      case PS_IN_MODEL:
        CHECK_START_ELE(command_line, PS_IN_COMMAND_LINE, ALL_CHARS);
        CHECK_START_ELE(max_correlation, PS_IN_MAX_CORRELATION, ALL_CHARS);
        DO_START_ELE(remove_correlated, PS_IN_REMOVE_CORRELATED, IGNORE);
        DO_START_ELE(strand_handling, PS_IN_STRAND_HANDLING, IGNORE);
        DO_START_ELE(translate_dna, PS_IN_TRANSLATE_DNA, IGNORE);
        CHECK_START_ELE(max_seq_evalue, PS_IN_MAX_SEQ_EVALUE, ALL_CHARS);
        DO_START_ELE(adj_hit_pvalue, PS_IN_ADJ_HIT_PVALUE, IGNORE);
        CHECK_START_ELE(max_hit_pvalue, PS_IN_MAX_HIT_PVALUE, ALL_CHARS);
        CHECK_START_ELE(max_weak_pvalue, PS_IN_MAX_WEAK_PVALUE, ALL_CHARS);
        CHECK_START_ELE(host, PS_IN_HOST, ALL_CHARS);
        CHECK_START_ELE(when, PS_IN_WHEN, ALL_CHARS);
        known = 0;
        break;
      case PS_IN_ALPHABET:
        DO_START_ELE(letter, PS_IN_LETTER, IGNORE);
        known = 0;
        break;
      case PS_IN_MOTIFS:
        DO_START_ELE(motif, PS_IN_MOTIF, IGNORE);
        DO_START_ELE(correlation, PS_IN_CORRELATION, IGNORE);
        DO_START_ELE(nos, PS_IN_NOS, IGNORE);
        known = 0;
        break;
      case PS_IN_NOS:
        DO_START_ELE(expect, PS_IN_EXPECT, IGNORE);
        known = 0;
        break;
      case PS_IN_SEQUENCES:
        DO_START_ELE(database, PS_IN_DATABASE, IGNORE);
        DO_START_ELE(sequence, PS_IN_SEQUENCE, IGNORE);
        known = 0;
        break;
      case PS_IN_SEQUENCE:
        DO_START_ELE(score, PS_IN_SCORE, IGNORE);
        DO_START_ELE(seg, PS_IN_SEG, IGNORE);
        known = 0;
        break;
      case PS_IN_SEG:
        CHECK_START_ELE(data, PS_IN_DATA, ALL_BUT_SPACE);
        DO_START_ELE(hit, PS_IN_HIT, IGNORE);
      default:
        known = 0;
    }
    if (!known) ps->udepth = 1;
  }
}

#define DO_END_ELE(_expected_,_state_,_transition_) \
  case _state_: \
    if (strcmp((char*)name, #_expected_) == 0) { \
      known = 1; \
      end_ele_ ## _expected_ (ps); \
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
  PARSER_STATE_T *ps = (PARSER_STATE_T*)ctx;
  int known;
  if (ps->state == PS_ERROR) return;
  if (ps->udepth) {
    ps->udepth -= 1; 
  } else {
    known = 0;
    switch (ps->state) {
      DO_END_ELE(mast, PS_IN_MAST, PS_END);
      CHECK_END_ELE(model, PS_IN_MODEL, PS_IN_MAST);
      DO_END_ELE(command_line, PS_IN_COMMAND_LINE, PS_IN_MODEL);
      DO_END_ELE(max_correlation, PS_IN_MAX_CORRELATION, PS_IN_MODEL);
      CHECK_END_ELE(remove_correlated, PS_IN_REMOVE_CORRELATED, PS_IN_MODEL);
      CHECK_END_ELE(strand_handling, PS_IN_STRAND_HANDLING, PS_IN_MODEL);
      CHECK_END_ELE(translate_dna, PS_IN_TRANSLATE_DNA, PS_IN_MODEL);
      DO_END_ELE(max_seq_evalue, PS_IN_MAX_SEQ_EVALUE, PS_IN_MODEL);
      CHECK_END_ELE(adj_hit_pvalue, PS_IN_ADJ_HIT_PVALUE, PS_IN_MODEL);
      DO_END_ELE(max_hit_pvalue, PS_IN_MAX_HIT_PVALUE, PS_IN_MODEL);
      DO_END_ELE(max_weak_pvalue, PS_IN_MAX_WEAK_PVALUE, PS_IN_MODEL);
      DO_END_ELE(host, PS_IN_HOST, PS_IN_MODEL);
      DO_END_ELE(when, PS_IN_WHEN, PS_IN_MODEL);
      CHECK_END_ELE(alphabet, PS_IN_ALPHABET, PS_IN_MAST);
      CHECK_END_ELE(letter, PS_IN_LETTER, PS_IN_ALPHABET);
      CHECK_END_ELE(motifs, PS_IN_MOTIFS, PS_IN_MAST);
      CHECK_END_ELE(motif, PS_IN_MOTIF, PS_IN_MOTIFS);
      CHECK_END_ELE(correlation, PS_IN_CORRELATION, PS_IN_MOTIFS);
      DO_END_ELE(nos, PS_IN_NOS, PS_IN_MOTIFS);
      CHECK_END_ELE(expect, PS_IN_EXPECT, PS_IN_NOS);
      DO_END_ELE(sequences, PS_IN_SEQUENCES, PS_IN_MAST);
      CHECK_END_ELE(database, PS_IN_DATABASE, PS_IN_SEQUENCES);
      DO_END_ELE(sequence, PS_IN_SEQUENCE, PS_IN_SEQUENCES);
      CHECK_END_ELE(score, PS_IN_SCORE, PS_IN_SEQUENCE);
      DO_END_ELE(seg, PS_IN_SEG, PS_IN_SEQUENCE);
      DO_END_ELE(data, PS_IN_DATA, PS_IN_SEG);
      CHECK_END_ELE(hit, PS_IN_HIT, PS_IN_SEG);
      CHECK_END_ELE(runtime, PS_IN_RUNTIME, PS_IN_MAST);

    }
    if (!known) {
      DEBUG_FMT(QUIET_VERBOSE,"Hit error at %s in state %d\n", (char*)name, ps->state);
      ps->state = PS_ERROR;
    }
  }
}

int parse_xml_file(const char *xml_file_name, char *txt_file_name) {
  PARSER_STATE_T my_state;
  xmlSAXHandler my_handler;

  DEBUG_MSG(NORMAL_VERBOSE, "Starting xml parse\n");

  my_state.txt_file_name = txt_file_name;
  memset(&my_handler, 0, sizeof(xmlSAXHandler));
  my_handler.startDocument = handle_start_doc;
  my_handler.endDocument = handle_end_doc;
  my_handler.characters = handle_characters;
  my_handler.startElement = handle_start_ele;
  my_handler.endElement = handle_end_ele;

  if (xmlSAXUserParseFile(&my_handler, &my_state, xml_file_name) < 0) {
    DEBUG_MSG(QUIET_VERBOSE, "Parser encountered structural errors\n");
    return -1;
  } else {
    if (my_state.state == PS_END) {
      DEBUG_MSG(NORMAL_VERBOSE, "Parser completed\n");
    } else {
      DEBUG_FMT(QUIET_VERBOSE, "Parser did not reach end state (%d)\n", my_state.state);
    }
    return 0;
  }
}


int main(int argc, char **argv) {
  int result;
  //Using SAX scan through the xml document
  //load the model parameters,
  //load the alphabet
  //load the motifs
  //if the mode is not separate strand scoring then just read down the file and generate outputs
  //otherwise... output the lowest of the two scores and
  //keep a ordered tree of scores which haven't been inserted yet and write out the data to file
  //and keep a pointer to it, 
  if (argc == 3) {
    result = parse_xml_file(argv[1], argv[2]);
  } else {
    fprintf(stdout, "mast2txt <mast.xml> <mast.txt>\n");
    result = 1;
  }

  return result;
}
/*}}}*/
