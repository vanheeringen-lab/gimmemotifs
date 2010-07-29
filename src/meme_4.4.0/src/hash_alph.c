/*
 * $Id: hash_alph.c 4278 2009-12-23 09:58:37Z james_johnson $
 * 
 * $Log$
 * Revision 1.2  2005/10/25 19:06:39  nadya
 * rm old macro for Header, all info is taken care of by Id and Log.
 *
 * Revision 1.1.1.1  2005/07/29 00:20:12  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

/***********************************************************************
*								       *
*	MEME++							       *
*	Copyright 1994, The Regents of the University of California    *
*	Author: Timothy L. Bailey				       *
*								       *
*	Modified: tlb 11-12-97;
***********************************************************************/
/* hash_alph */

#include "user.h"
#include "macros.h"
#include "hash_alph.h"

/*
	Routines for hashing on the alphabet of sequences.
*/

/**********************************************************************/
/* 
	setup_hash_alph

	Sets up the hash table for the given alphabet.  There are three
	hash tables:
		x	alphabet
		--	--------
		0:	DNAB
		1:	PROTEINB
		2:	other (including DNA0 and PROTEIN0)
	This function can be called once for each of these alphabets,
	and then the current alphabet can be selected using setalph(x).

	Unused letters hash to -1.

	unhash(alength) will give 'X' (only works for alphabet set up last).

	If the alphabet is DNA0 or PROTEIN0, all ambiguous characters in
	the BLAST alphabet are hashed to 'X'.
*/
/**********************************************************************/
extern int setup_hash_alph(
  char *alphabet			/* the alphabet to set up hashing for */
)
{
  int i, alength;
  char c;

  /* 
   get length of alphabet 
  */
  alength = strlen(alphabet);
  if (alength > MAXALPH) {
    fprintf(stderr, "Alphabet too long (> %d).\n", MAXALPH);
    exit(1);
  }

  /*
    determine if current alphabet DNAB or PROTEINB or other
  */
  if (!strcmp(alphabet, DNAB)) {
    setalph(0);					/* current alphabet DNAB */
  } else if (!strcmp(alphabet, PROTEINB)) {
    setalph(1);					/* current alphabet PROTEINB */
  } else {
    setalph(2);					/* current alphabet other */
  }

  /* 
    flag unused letters 
  */
  for (i=0; i<MAXASCII; i++) hash(i) = -1;

  /* 
    set up the hashing and unhashing indices
  */
  for (i = 0; (c = alphabet[i]); i++) {
    c = islower((int)c) ? toupper((int)c) : c;	/* convert to uppercase */
    hash(c) = i;
    unhash(i) = c;
  }

  /* 
    MEME: convert ambiguous characters to X
  */
  if (!strcmp(alphabet, DNA0)) { 		/* DNA */
    hash('X') = alength;
    unhash(alength) = 'X';
    for (i=0; (c=DNAB[i]); i++) if (!strchr(DNA0, c)) 
      hash(c) = alength;
  } else if (!strcmp(alphabet, PROTEIN0)) { 	/* PROTEIN */
    hash('X') = alength;
    unhash(alength) = 'X';
    for (i=0; (c=PROTEINB[i]); i++) if (!strchr(PROTEIN0, c)) 
      hash(c) = alength;
  }

  return alength;
} /* setup_hash_alph */

/**********************************************************************/
/*
	r2seq

	Convert integer-coded sequence to ascii.
*/
/**********************************************************************/
extern void r2seq(
  char *seq,
  char *res,
  int len
)
{
  int i;

  for (i=0; i<len; i++) {
    seq[i] = unhash(res[i]);
  }
  seq[i] = '\0';
}

/***********************************************************************/
/*
	get_blast_alphabet

	Find out which BLAST alphabet it is a subset of, extend
	it to the BLAST alphabet and create a permutation matrix
	that says which letter(s) in the old alphabet should be averaged
	together to create each position in the new alphabet.  Rows
	of permutation matrix are terminated with -1 (like NULL in string).

	Returns the new alphabet.
	
	Updates p.
*/
/***********************************************************************/
extern char *get_blast_alphabet(
  char *old_alph, 		/* old alphabet IN (but converted to uppercase) */
  int *p[MAXASCII]		/* permutation and substitution matrix OUT */
)
{
  int i, j;
  int old_alen = strlen(old_alph);	/* length of old alphabet */
  char *new_alph;			/* new (complete) alphabet */
  int new_alen;				/* length of new alphabet */
  char *to;				/* list of substitutions */
  int to_len;				/* length of to */
  char **subst;				/* substitution list */

  /* convert alphabet to uppercase */
  for (i=0; i<old_alen; i++) 
    if (islower((int)old_alph[i])) old_alph[i] = toupper((int)old_alph[i]);

  /* 
    determine what type of alphabet we have 
  */
  if (strspn(old_alph, DNAB) == old_alen) {
    setalph(0);			/* set alphabet hash function */
    new_alph = DNAB; 		/* BLAST DNA alphabet */
    subst = dna_subst;		/* list of substitutions */
  } else if (strspn(old_alph, PROTEINB) == old_alen) {
    setalph(1);			/* set alphabet hash function */
    new_alph = PROTEINB;	/* BLAST PROTEIN alphabet */
    subst = prot_subst;		/* list of substitutions */
  } else {
    fprintf(stderr, "Don't recognize the motif alphabet: %s\n", old_alph);
    exit(1);
  }

  /* 
    create mapping from new position to old positions 
  */
  new_alen = strlen(new_alph);		/* length of new alphabet */
  for (i=0; i < new_alen; i++) { 	/* init mapping matrix */
    p[i] = NULL; 			/* create row of substitutions */
    Resize(p[i], 2, int); 		/* make row length 2 */
    p[i][0] = -1;			/* flag end of substitution list */
  }
  for (i=0; i < new_alen; i++) {	/* new position */
    char c = new_alph[i];		/* new letter */
    char *o;				/* pointer to c in old alphabet */ 
    if ((o = strchr(old_alph, c)) != NULL) {	/* letter is in old alphabet */
      p[i][0] = (int)(o-old_alph);	/* put in list */
      p[i][1] = -1;			/* flag end of list */ 
    } else {				/* letter not in old alphabet */
      to = subst[i];			/* letters to substitute */
      to_len = strlen(to);		/* number of substitutions */
      Resize(p[i], to_len+1, int);	/* make list of substitions */
      for (j=0; j<to_len; j++) {
        if ((o = strchr(old_alph, to[j])) != NULL) {	/* letter found */
          p[i][j] = (int) (o - old_alph);	/* put in list */
        } else {			/* required letter missing */
          char *a = (subst==dna_subst ? "DNA" : "protein");
          fprintf(stderr, 
            "The motif alphabet %s appears to be a %s alphabet\n", old_alph, a);
 	  fprintf(stderr, "but is missing the required letter `%c'.\n", to[j]);
          exit(1);
        }
      }
      p[i][j] = -1;			/* flag end of list */ 
    }
  }					/* new position */

  return new_alph;
} /* get_blast_alphabet */

/**********************************************************************/
/*
	setup_hash_dna2prot

	Set up hash table for hashing from DNAB letters to PROTEINB position.

	Routine setup_hash must be called first to setup hashing functions
	dnabhash and protbhash.
	
*/
/**********************************************************************/
extern void setup_hash_dnab2protb()
{
  int i, j;
  int i0, i1, i2;		/* index in DNAB alphabet of codon */
  int j0, j1, j2;		/* positions in substitution lists */
  char s0, s1, s2;		/* substitution characters (DNA0 alphabet) */
  int subs[MAXASCII];		/* L (PROTEIN0) in map -> subs[protbhash(L)]=1*/
  char substr[MAXALPH+1];	/* string of substitution letters for codon */
  char *protalph = PROTEINB;	/* PROTEINB alphabet */
  int protalen = strlen(PROTEINB);/* length of PROTEINB alphabet */
  char *dnalph = DNAB;		/* DNAB alphabet */
  char aa;			/* amino acid letter */

  /* set length of DNAB alphabet for use by macro dnab2protb */
  dnablen = strlen(DNAB);
  
  /* create the dnab2protb index table */
  dnab2protb_index = NULL;
  Resize(dnab2protb_index, dnablen*dnablen*dnablen, int);

  /* 
    For each possible DNAB codon, find all the translations of it
    (more than one if it contains ambiguous DNAB characters)
    and see if if matches one of the PROTEINB characters.  Otherwise,
    map it to "X".  Store the position in PROTEINB alphabet in the index table.
  */
  for (i0=0; i0<dnablen; i0++) {			/* codon position 0 */	
    for (i1=0; i1<dnablen; i1++) {			/* codon position 1 */
      for (i2=0; i2<dnablen; i2++) {			/* codon position 2 */
        for (i=0; i<protalen; i++) subs[i] = 0;		/* clear sub list */
	for (j0=0; (s0=dna_subst[i0][j0]); j0++) {	/* position 0 subst. */
	  for (j1=0; (s1=dna_subst[i1][j1]); j1++) {	/* position 1 subst. */
	    for (j2=0; (s2=dna_subst[i2][j2]); j2++) {	/* position 2 subst. */
	      aa = dna0_to_prot0(s0, s1, s2);     	/* translated aa */
              subs[protbhash(aa)] = 1;			/* mark substitution */
            } /* subst 2 */
          } /* subst 1 */
        } /* subst 0 */

        /* 
          All substitutions of current codon (i0, i1, i2) have now been
	  marked in the subs list.  Convert the list to a string and
          find the PROTEINB character whose substitution list matches it.
          Otherwise, set the match to "X".
        */
        j = 0;				/* number of substitution letters */
        for (i=0; i<protalen; i++) if (subs[i]) substr[j++] = protalph[i]; 
        if (j==1) {			/* unambiguous codon */
          aa = substr[0];
        } else if (j==protalen) {	/* full substitution list */
          aa = 'X';
        } else {			/* ambiguous codon */
          substr[j] = '\0';		/* NULL at end of string */
	  /* 
            Search for substitution list that matches this codons subst list.
           */
          aa = 'X';			/* default if no match found */
          for (i=0; i<protalen; i++) {	/* position in PROTEINB alphabet */
            if (!strcmp(substr, prot_subst[i])) {	/* match */
              aa = protalph[i];
              break;
            }
          } /* position in PROTEINB alphabet */
         } /* ambiguous codon */

         /* set the lookup table entry for this codon */
         hash_dnab2protb(dnalph[i0],dnalph[i1],dnalph[i2]) = protbhash(aa);
      } /* codon position 2 */
    } /* codon position 1 */
  } /* codon position 0 */
} /* setup_hash_dnab2protb */

/**********************************************************************/
/*
        invcomp_dna
 
        Convert a DNA sequence in place to its reverse complement.
*/
/**********************************************************************/
extern void invcomp_dna(
  char *sequence,                       /* DNA sequence */
  long length                            /* length of sequence */
)
{
  char *sl = sequence;                  /* left end of sequence */
  char *sr = sequence + length - 1;     /* right end of sequence */
  for (; sl<=sr; sl++, sr--) {
    char tmp = comp_dna(*sl);
    *sl = comp_dna(*sr);
    *sr = tmp;
  }
} /* invcomp_dna */

/**********************************************************************/
/*
	dhash_it

	Hash sequence to index in two-letter hashed alphabet.
	If translating DNA, hash two codons to two-letter index.
	Compute the letter frequencies in the sequence.
	If translating DNA, report frequencies of letters translated
	in reading frame with fewest stop codons.

*/
/**********************************************************************/
extern int *dhash_it(
  BOOLEAN xlate_dna,			/* database is DNA and motifs protein */
  int alen,				/* length of alphabet */
  char *sequence,			/* sequence of sample */
  long length 				/* length of sequence */
)
{
  long i;
  int *hash_seq = NULL;			/* hashed sequence */
  int *h;				/* pointer in hash_seq */
  char *s;				/* pointer in sequence */
  long len;				/* length of sequence */
  int inc;				/* distance to next letter in hashed */

  /*
    Hash sequence from letters to positions in alphabet.
    Three adjacent letters (codons) are hashed if translating DNA. 
  */
  Resize(hash_seq, length+3, int);		/* leave room for padding */
  len = xlate_dna ? length - 2 : length;	/* last full letter or codon */
  for (i=0,s=sequence,h=hash_seq; i<len; i++,s++,h++) 
    *h = chash(xlate_dna, FALSE, s);
  for ( ; i<len+3; i++,h++) *h = alen;		/* pad with alen */

  /*
    Hash sequence to "double letter" logodds alphabet.
  */
  inc = xlate_dna ? 3 : 1;			/* distance to next letter */
  for (i=0, h=hash_seq; i<len; i++,h++) *h = dhash(*h, *(h+inc), alen);

  return hash_seq;
} /* dhash_it */

