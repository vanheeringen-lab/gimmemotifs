/*
 * $Id: ureadseq.c 4080 2009-10-19 20:18:09Z cegrant $
 * 
 * $Log$
 * Revision 1.2  2006/03/08 20:50:11  nadya
 * merge chamges from v3_5_2 branch
 *
 * Revision 1.1.1.1.4.1  2006/01/26 08:34:26  tbailey
 * Renamed local function getline() to getline1() to avoid conflict
 * with system function defined in stdio.h
 *
 * Revision 1.1.1.1  2005/07/29 17:19:22  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

/* File: ureadseq.c
 *
 * Reads and writes nucleic/protein sequence in various
 * formats. Data files may have multiple sequences.
 *
 * Copyright 1990 by d.g.gilbert
 * biology dept., indiana university, bloomington, in 47405
 * e-mail: gilbertd@bio.indiana.edu
 *
 * This program may be freely copied and used by anyone.
 * Developers are encourged to incorporate parts in their
 * programs, rather than devise their own private sequence
 * format.
 *
 * This should compile and run with any ANSI C compiler.
 *
 */


#define UREADSEQ_G
#include "ureadseq.h"

/* strlcpy is missing from some LINUX */
#if defined(Linux)
static size_t strlcpy(char *dst, const char *src, size_t dstsize)
{
  int i;
  for (i=0; src[i] != '\0'; i++) {
    if (i<dstsize) dst[i] = src[i];
  }
  if (i<dstsize) dst[i] = '\0'; else dst[dstsize-1] = '\0';
  return(i); 
}
#endif

int Strcasecmp(const char *a, const char *b)  /* from Nlm_StrICmp */
{
  int diff, done;
  if (a == b)  return 0;
  done = 0;
  while (! done) {
    diff = to_upper(*a) - to_upper(*b);
    if (diff) return diff;
    if (*a == '\0') done = 1;
    else { a++; b++; }
    }
  return 0;
}

int Strncasecmp(const char *a, const char *b, long maxn) /* from Nlm_StrNICmp */
{
  int diff, done;
  if (a == b)  return 0;
  done = 0;
  while (! done) {
    diff = to_upper(*a) - to_upper(*b);
    if (diff) return diff;
    if (*a == '\0') done = 1;
    else {
      a++; b++; maxn--;
      if (! maxn) done = 1;
      }
    }
  return 0;
}





#ifndef Local
# define Local      static    /* local functions */
#endif

#define kStartLength  500

const char *aminos      = "ABCDEFGHIKLMNPQRSTVWXYZ*";
const char *primenuc    = "ACGTU";
const char *protonly    = "EFIPQZ";

const char kNocountsymbols[5]  = "_.-?";
const char stdsymbols[6]  = "_.-*?";
const char allsymbols[32] = "_.-*?<>{}[]()!@#$%^&=+;:'/|`~\"\\";
static const char *seqsymbols   = allsymbols;

const char nummask[11]   = "0123456789";
const char nonummask[11] = "~!@#$%^&*(";

/*
    use general form of isseqchar -- all chars + symbols.
    no formats except nbrf (?) use symbols in data area as
    anything other than sequence chars.
*/



                          /* Local variables for readSeq: */
struct ReadSeqVars {
  short choice, err, nseq;
  long  seqlen, maxseq, seqlencount;
  short topnseq;
  long  topseqlen;
  const char *fname;
  char *seq, *seqid, matchchar;
  boolean allDone, done, filestart, addit;
  FILE  *f;
  long  linestart;
  char  s[MAXLINE], *sp;

  int  (*isseqchar)(int c);	/* tlb 3/1/96 */
  /*int (*isseqchar)();*/
  /* int  (*isseqchar)(int c);  << sgi cc hates (int c) */
};



int isSeqChar(int c)
{
  return (isalpha((int)c) || strchr(seqsymbols,c));
}

int isSeqNumChar(int c)
{
  return (isalnum((int)c) || strchr(seqsymbols,c));
}


int isAnyChar(int c)
{
  return isascii(c); /* wrap in case isascii is macro */
}

Local void readline(FILE *f, char *s, long *linestart)
{
  char  *cp;

  *linestart= ftell(f);
  if (NULL == fgets(s, MAXLINE, f))
    *s = 0;
  else {
    cp = strchr(s, '\n');
    if (cp != NULL) *cp = 0;
    }
}

Local void getline1(struct ReadSeqVars *V)
{
  readline(V->f, V->s, &V->linestart);
}

Local void ungetline(struct ReadSeqVars *V)
{
  fseek(V->f, V->linestart, 0);
}


Local void addseq(char *s, struct ReadSeqVars *V)
{
  char  *ptr;

  if (V->addit) while (*s != 0) {
    if ((V->isseqchar)(*s)) {
      if (V->seqlen >= V->maxseq) {
        V->maxseq += kStartLength;
        ptr = (char*) realloc(V->seq, V->maxseq+1);
        if (ptr==NULL) {
          V->err = eMemFull;
          return;
          }
        else V->seq = ptr;
        }
      V->seq[(V->seqlen)++] = *s;
      }
    s++;
    }
}

Local void countseq(char *s, struct ReadSeqVars *V)
 /* this must count all valid seq chars, for some formats (paup-sequential) even
    if we are skipping seq... */
{
  while (*s != 0) {
    if ((V->isseqchar)(*s)) {
      (V->seqlencount)++;
      }
    s++;
    }
}


Local void addinfo(char *s, struct ReadSeqVars *V)
{
  char *si = (char *) malloc((strlen(s) + 40) * sizeof(char));
  boolean saveadd;

  while (*s == ' ') s++;
  sprintf(si, " %d)  %s\n", V->nseq, s);

  saveadd = V->addit;
  V->addit = true;
  V->isseqchar = isAnyChar;
  addseq( si, V);
  free(si);
  V->addit = saveadd;
  V->isseqchar = isSeqChar;
}




Local void readLoop(short margin, boolean addfirst,
            boolean (*endTest)(boolean *addend, boolean *ungetend, struct ReadSeqVars *V),
            struct ReadSeqVars *V)
{
  boolean addend = false;
  boolean ungetend = false;

  V->nseq++;
  if (V->choice == kListSequences) V->addit = false;
  else V->addit = (V->nseq == V->choice);
  if (V->addit) V->seqlen = 0;

  if (addfirst) addseq(V->s, V);
  do {
    getline1(V);
    V->done = feof(V->f);
    V->done |= (*endTest)( &addend, &ungetend, V);
    if (V->addit && (addend || !V->done) && ((int) strlen(V->s) > margin)) {
      addseq( (V->s)+margin, V);
    }
  } while (!V->done);

  if (V->choice == kListSequences) addinfo(V->seqid, V);
  else {
    V->allDone = (V->nseq >= V->choice);
    if (V->allDone && ungetend) ungetline(V);
    }
}



Local boolean endIG( boolean *addend, boolean *ungetend, struct ReadSeqVars *V)
{
  *addend = true; /* 1 or 2 occur in line w/ bases */
  *ungetend= false;
  return((strchr(V->s,'1')!=NULL) || (strchr(V->s,'2')!=NULL));
}

Local void readIG(struct ReadSeqVars *V)
{
/* 18Aug92: new IG format -- ^L between sequences in place of ";" */
  char  *si;

  while (!V->allDone) {
    do {
      getline1(V);
      for (si= V->s; *si != 0 && *si < ' '; si++) *si= ' '; /* drop controls */
      if (*si == 0) *V->s= 0; /* chop line to empty */
    } while (! (feof(V->f) || ((*V->s != 0) && (*V->s != ';') ) ));
    if (feof(V->f))
      V->allDone = true;
    else {
      strcpy(V->seqid, V->s);
      readLoop(0, false, endIG, V);
      }
  }
}



Local boolean endStrider( boolean *addend, boolean *ungetend, struct ReadSeqVars *V)
{
  *addend = false;
  *ungetend= false;
  return (strstr( V->s, "//") != NULL);
}

Local void readStrider(struct ReadSeqVars *V)
{ /* ? only 1 seq/file ? */

  while (!V->allDone) {
    getline1(V);
    if (strstr(V->s,"; DNA sequence  ") == V->s)
      strcpy(V->seqid, (V->s)+16);
    else
      strcpy(V->seqid, (V->s)+1);
    while ((!feof(V->f)) && (*V->s == ';')) {
      getline1(V);
      }
    if (feof(V->f)) V->allDone = true;
    else readLoop(0, true, endStrider, V);
  }
}


Local boolean endPIR( boolean *addend, boolean *ungetend, struct ReadSeqVars *V)
{
  *addend = false;
  *ungetend= (strstr(V->s,"ENTRY") == V->s);
  return ((strstr(V->s,"///") != NULL) || *ungetend);
}

Local void readPIR(struct ReadSeqVars *V)
{ /*PIR -- many seqs/file */

  while (!V->allDone) {
    while (! (feof(V->f) || strstr(V->s,"ENTRY")  || strstr(V->s,"SEQUENCE")) )
      getline1(V);
    strcpy(V->seqid, (V->s)+16);
    while (! (feof(V->f) || strstr(V->s,"SEQUENCE") == V->s))
      getline1(V);
    readLoop(0, false, endPIR, V);

    if (!V->allDone) {
     while (! (feof(V->f) || ((*V->s != 0)
       && (strstr( V->s,"ENTRY") == V->s))))
        getline1(V);
      }
    if (feof(V->f)) V->allDone = true;
  }
}


Local boolean endGB( boolean *addend, boolean *ungetend, struct ReadSeqVars *V)
{
  *addend = false;
  *ungetend= (strstr(V->s,"LOCUS") == V->s);
  return ((strstr(V->s,"//") != NULL) || *ungetend);
}

Local void readGenBank(struct ReadSeqVars *V)
{ /*GenBank -- many seqs/file */

  while (!V->allDone) {
    strcpy(V->seqid, (V->s)+12);
    while (! (feof(V->f) || strstr(V->s,"ORIGIN") == V->s))
      getline1(V);
    readLoop(0, false, endGB, V);

    if (!V->allDone) {
     while (! (feof(V->f) || ((*V->s != 0)
       && (strstr( V->s,"LOCUS") == V->s))))
        getline1(V);
      }
    if (feof(V->f)) V->allDone = true;
  }
}


Local boolean endNBRF( boolean *addend, boolean *ungetend, struct ReadSeqVars *V)
{
  char  *a;

  if ((a = strchr(V->s, '*')) != NULL) { /* end of 1st seq */
    /* "*" can be valid base symbol, drop it here */
    *a = 0;
    *addend = true;
    *ungetend= false;
    return(true);
    }
  else if (*V->s == '>') { /* start of next seq */
    *addend = false;
    *ungetend= true;
    return(true);
    }
  else
    return(false);
}

Local void readNBRF(struct ReadSeqVars *V)
{
  while (!V->allDone) {
    strcpy(V->seqid, (V->s)+4);
    getline1(V);   /*skip title-junk line*/
    readLoop(0, false, endNBRF, V);
    if (!V->allDone) {
     while (!(feof(V->f) || (*V->s != 0 && *V->s == '>')))
        getline1(V);
      }
    if (feof(V->f)) V->allDone = true;
  }
}



Local boolean endPearson( boolean *addend, boolean *ungetend, struct ReadSeqVars *V)
{
  *addend = false;
  *ungetend= true;
  return(*V->s == '>');
}

Local void readPearson(struct ReadSeqVars *V)
{
  while (!V->allDone) {
    strlcpy(V->seqid, (V->s)+1, MAXLINE);
    readLoop(0, false, endPearson, V);
    if (!V->allDone) {
     while (!(feof(V->f) || ((*V->s != 0) && (*V->s == '>'))))
        getline1(V);
      }
    if (feof(V->f)) V->allDone = true;
  }
}



Local boolean endEMBL( boolean *addend, boolean *ungetend, struct ReadSeqVars *V)
{
  *addend = false;
  *ungetend= (strstr(V->s,"ID   ") == V->s);
  return ((strstr(V->s,"//") != NULL) || *ungetend);
}

Local void readEMBL(struct ReadSeqVars *V)
{
  while (!V->allDone) {
    strcpy(V->seqid, (V->s)+5);
    do {
      getline1(V);
    } while (!(feof(V->f) | (strstr(V->s,"SQ   ") == V->s)));

    readLoop(0, false, endEMBL, V);
    if (!V->allDone) {
      while (!(feof(V->f) |
         ((*V->s != '\0') & (strstr(V->s,"ID   ") == V->s))))
      getline1(V);
    }
    if (feof(V->f)) V->allDone = true;
  }
}



Local boolean endZuker( boolean *addend, boolean *ungetend, struct ReadSeqVars *V)
{
  *addend = false;
  *ungetend= true;
  return( *V->s == '(' );
}

Local void readZuker(struct ReadSeqVars *V)
{
  /*! 1st string is Zuker's Fortran format */

  while (!V->allDone) {
    getline1(V);  /*s == "seqLen seqid string..."*/
    strcpy(V->seqid, (V->s)+6);
    readLoop(0, false, endZuker, V);
    if (!V->allDone) {
      while (!(feof(V->f) |
        ((*V->s != '\0') & (*V->s == '('))))
          getline1(V);
      }
    if (feof(V->f)) V->allDone = true;
  }
}



Local boolean endFitch( boolean *addend, boolean *ungetend, struct ReadSeqVars *V)
{
  /* this is a somewhat shaky end,
    1st char of line is non-blank for seq. title
  */
  *addend = false;
  *ungetend= true;
  return( *V->s != ' ' );
}

Local void readFitch(struct ReadSeqVars *V)
{
  boolean first;

  first = true;
  while (!V->allDone) {
    if (!first) strcpy(V->seqid, V->s);
    readLoop(0, first, endFitch, V);
    if (feof(V->f)) V->allDone = true;
    first = false;
    }
}


Local void readPlain(struct ReadSeqVars *V)
{
  V->nseq++;
  V->addit = (V->choice > 0);
  if (V->addit) V->seqlen = 0;
  addseq(V->seqid, V);   /*from above..*/
  if (V->fname!=NULL) sprintf(V->seqid, "%s  [Unknown form]", V->fname);
  else sprintf(V->seqid, "  [Unknown form]");
  do {
    addseq(V->s, V);
    V->done = feof(V->f);
    getline1(V);
  } while (!V->done);
  if (V->choice == kListSequences) addinfo(V->seqid, V);
  V->allDone = true;
}


Local void readUWGCG(struct ReadSeqVars *V)
{
/*
10nov91: Reading GCG files casued duplication of last line when
         EOF followed that line !!!
    fix: getline1 now sets *V->s = 0
*/
  char  *si;

  V->nseq++;
  V->addit = (V->choice > 0);
  if (V->addit) V->seqlen = 0;
  strcpy(V->seqid, V->s);
  /*writeseq: "    %s  Length: %d  (today)  Check: %d  ..\n" */
  /*drop above or ".." from id*/
  if ( (si = strstr(V->seqid,"  Length: ")) ) *si = 0;
  else if ( (si = strstr(V->seqid,"..")) ) *si = 0;
  do {
    V->done = feof(V->f);
    getline1(V);
    if (!V->done) addseq((V->s), V);
  } while (!V->done);
  if (V->choice == kListSequences) addinfo(V->seqid, V);
  V->allDone = true;
}


Local void readOlsen(struct ReadSeqVars *V)
{ /* G. Olsen /print output from multiple sequence editor */

  char    *si, *sj, *sk, *sm=NULL, sid[40], snum[20];
  boolean indata = false;
  int snumlen = 0;

  V->addit = (V->choice > 0);
  if (V->addit) V->seqlen = 0;
  rewind(V->f); V->nseq= 0;
  do {
    getline1(V);
    V->done = feof(V->f);

    if (V->done && !(*V->s)) break;
    else if (indata) {
      if ( (si= strstr(V->s, sid))
        /* && (strstr(V->s, snum) == si - snumlen - 1) ) { */
        && (sm= strstr(V->s, snum)) && (sm < si - snumlen) ) {

        /* Spaces are valid alignment data !! */
/* 17Oct91: Error, the left margin is 21 not 22! */
/* dropped some nucs up to now -- my example file was right shifted ! */
/* variable right id margin, drop id-2 spaces at end */
/*
  VMS CC COMPILER (VAXC031) mess up:
  -- Index of 21 is chopping 1st nuc on VMS systems Only!
  Byte-for-byte same ame rnasep.olsen sequence file !
*/

        /* si = (V->s)+21; < was this before VMS CC wasted my time */
        si += 10;  /* use strstr index plus offset to outfox VMS CC bug */

        if ( (sk = strstr(si, sid)) ) *(sk-2) = 0;
        for (sk = si; *sk != 0; sk++) {
           if (*sk == ' ') *sk = '.';
           /* 18aug92: !! some olsen masks are NUMBERS !! which addseq eats */
           else if (isdigit((int)*sk)) *sk= nonummask[*sk - '0'];
           }

        addseq(si, V);
        }
      }

    else if ( (sk = strstr(V->s, "): ")) ) {  /* seq info header line */
  /* 18aug92: correct for diff seqs w/ same name -- use number, e.g. */
  /*   3 (Agr.tume):  agrobacterium.prna  18-JUN-1987 16:12 */
  /* 328 (Agr.tume):  agrobacterium.prna XYZ  19-DEC-1992   */
      (V->nseq)++;
      si = 1 + strchr(V->s,'(');
      *sk = ' ';
      if (V->choice == kListSequences) addinfo( si, V);
      else if (V->nseq == V->choice) {
        strcpy(V->seqid, si);
        sj = strchr(V->seqid, ':');
        while (*(--sj) == ' ') ;
        while (--sj != V->seqid) { if (*sj == ' ') *sj = '_'; }

        *sk = 0;
        while (*(--sk) == ' ') *sk = 0;
        strcpy(sid, si);

        si= V->s;
        while ((*si <= ' ') && (*si != 0)) si++;
        snumlen=0;
        while (si[snumlen] > ' ' && snumlen<20)
         { snum[snumlen]= si[snumlen]; snumlen++; }
        snum[snumlen]= 0;
        }

      }

    else if (strstr(V->s,"identity:   Data:")) {
      indata = true;
      if (V->choice == kListSequences) V->done = true;
      }

  } while (!V->done);

  V->allDone = true;
} /*readOlsen*/


Local void readMSF(struct ReadSeqVars *V)
{ /* gcg's MSF, mult. sequence format, interleaved ! */

  char    *si, *sj, sid[128];
  boolean indata = false;
  int     iline= 0;

  V->addit = (V->choice > 0);
  if (V->addit) V->seqlen = 0;
  rewind(V->f); V->nseq= 0;
  do {
    getline1(V);
    V->done = feof(V->f);

    if (V->done && !(*V->s)) break;
    else if (indata) {
      /*somename  ...gpvedai .......t.. aaigr..vad tvgtgptnse aipaltaaet */
      /*       E  gvenae.kgv tentna.tad fvaqpvylpe .nqt...... kv.affynrs */

      si= V->s;
      skipwhitespace(si);
      /* for (sj= si; isalnum((int)*sj); sj++) ; bug -- cdelwiche uses "-", "_" and others in names*/
      for (sj= si; *sj > ' '; sj++) ;
      *sj= 0;
      if ( *si ) {
        if ( (0==strcmp(si, sid)) ) {
          addseq(sj+1, V);
          }
        iline++;
        }
      }

    else if (NULL != (si = strstr(V->s, "Name: "))) {  /* seq info header line */
      /* Name: somename      Len:   100  Check: 7009  Weight:  1.00 */

      (V->nseq)++;
      si += 6;
      if (V->choice == kListSequences) addinfo( si, V);
      else if (V->nseq == V->choice) {
        strcpy(V->seqid, si);
        si = V->seqid;
        skipwhitespace(si);
        /* for (sj= si; isalnum((int)*sj); sj++) ; -- bug */
        for (sj= si; *sj > ' '; sj++) ;
        *sj= 0;
        strcpy(sid, si);
        }
      }

    else if ( strstr(V->s,"//") /*== V->s*/ )  {
      indata = true;
      iline= 0;
      if (V->choice == kListSequences) V->done = true;
      }

  } while (!V->done);


  V->allDone = true;
} /*readMSF*/



Local void readPAUPinterleaved(struct ReadSeqVars *V)
{ /* PAUP mult. sequence format, interleaved or sequential! */

  char    *si, *sj, *send, sid[40], sid1[40], saveseq[255];
  boolean first = true, indata = false, domatch;
  int     iline= 0, ifmc, saveseqlen=0;

#define fixmatchchar(s) { \
  for (ifmc=0; ifmc<saveseqlen; ifmc++) \
    if (s[ifmc] == V->matchchar) s[ifmc]= saveseq[ifmc]; }

  V->addit = (V->choice > 0);
  V->seqlencount = 0;
  if (V->addit) V->seqlen = 0;
  /* rewind(V->f); V->nseq= 0;  << do in caller !*/
  indata= true; /* call here after we find "matrix" */
  domatch= (V->matchchar > 0);

  do {
    getline1(V);
    V->done = feof(V->f);

    if (V->done && !(*V->s)) break;
    else if (indata) {
      /* [         1                    1                    1         ]*/
      /* human     aagcttcaccggcgcagtca ttctcataatcgcccacggR cttacatcct*/
      /* chimp     ................a.t. .c.................a ..........*/
      /* !! need to correct for V->matchchar */
      si= V->s;
      skipwhitespace(si);
      if (strchr(si,';')) indata= false;

      if (isalnum((int)*si))  {
        /* valid data line starts w/ a left-justified seq name in columns [0..8] */
        if (first) {
          (V->nseq)++;
          if (V->nseq >= V->topnseq) first= false;
          for (sj = si; isalnum((int)*sj); sj++) ;
          send= sj;
          skipwhitespace(sj);
          if (V->choice == kListSequences) {
            *send= 0;
            addinfo( si, V);
            }
          else if (V->nseq == V->choice) {
            if (domatch) {
              if (V->nseq == 1) { strcpy( saveseq, sj); saveseqlen= strlen(saveseq); }
              else fixmatchchar( sj);
              }
            addseq(sj, V);
            *send= 0;
            strcpy(V->seqid, si);
            strcpy(sid, si);
            if (V->nseq == 1) strcpy(sid1, sid);
            }
          }

        else if ( (strstr(si, sid) == si) ){
          while (isalnum((int)*si)) si++;
          skipwhitespace(si);
          if (domatch) {
            if (V->nseq == 1) { strcpy( saveseq, si); saveseqlen= strlen(saveseq); }
            else fixmatchchar( si);
            }
          addseq(si, V);
          }

        else if (domatch && (strstr(si, sid1) == si)) {
          strcpy( saveseq, si);
          saveseqlen= strlen(saveseq);
          }

        iline++;
        }
      }

    else if ( strstr(V->s,"matrix") )  {
      indata = true;
      iline= 0;
      if (V->choice == kListSequences) V->done = true;
      }

  } while (!V->done);

  V->allDone = true;
} /*readPAUPinterleaved*/



Local void readPAUPsequential(struct ReadSeqVars *V)
{ /* PAUP mult. sequence format, interleaved or sequential! */
  char    *si, *sj;
  boolean atname = true, indata = false;

  V->addit = (V->choice > 0);
  if (V->addit) V->seqlen = 0;
  V->seqlencount = 0;
  /* rewind(V->f); V->nseq= 0;  << do in caller !*/
  indata= true; /* call here after we find "matrix" */
  do {
    getline1(V);
    V->done = feof(V->f);

    if (V->done && !(*V->s)) break;
    else if (indata) {
      /* [         1                    1                    1         ]*/
      /* human     aagcttcaccggcgcagtca ttctcataatcgcccacggR cttacatcct*/
      /*           aagcttcaccggcgcagtca ttctcataatcgcccacggR cttacatcct*/
      /* chimp     ................a.t. .c.................a ..........*/
      /*           ................a.t. .c.................a ..........*/

      si= V->s;
      skipwhitespace(si);
      if (strchr(si,';')) indata= false;
      if (isalnum((int)*si))  {
        /* valid data line starts w/ a left-justified seq name in columns [0..8] */
        if (atname) {
          (V->nseq)++;
          V->seqlencount = 0;
          atname= false;
          sj= si+1;
          while (isalnum((int)*sj)) sj++;
          if (V->choice == kListSequences) {
            /* !! we must count bases to know when topseqlen is reached ! */
            countseq(sj, V);
            if (V->seqlencount >= V->topseqlen) atname= true;
            *sj= 0;
            addinfo( si, V);
            }
          else if (V->nseq == V->choice) {
            addseq(sj, V);
            V->seqlencount= V->seqlen;
            if (V->seqlencount >= V->topseqlen) atname= true;
            *sj= 0;
            strcpy(V->seqid, si);
            }
          else {
            countseq(sj, V);
            if (V->seqlencount >= V->topseqlen) atname= true;
            }
          }

        else if (V->nseq == V->choice) {
          addseq(V->s, V);
          V->seqlencount= V->seqlen;
          if (V->seqlencount >= V->topseqlen) atname= true;
          }
        else {
          countseq(V->s, V);
          if (V->seqlencount >= V->topseqlen) atname= true;
          }
        }
      }

    else if ( strstr(V->s,"matrix") )  {
      indata = true;
      atname= true;
      if (V->choice == kListSequences) V->done = true;
      }

  } while (!V->done);

  V->allDone = true;
} /*readPAUPsequential*/


Local void readPhylipInterleaved(struct ReadSeqVars *V)
{
  char    *si, *sj;
  boolean first = true;
  int     iline= 0;

  V->addit = (V->choice > 0);
  if (V->addit) V->seqlen = 0;
  V->seqlencount = 0;
  /* sscanf( V->s, "%d%d", &V->topnseq, &V->topseqlen); << topnseq == 0 !!! bad scan !! */
  si= V->s;
  skipwhitespace(si);
  V->topnseq= atoi(si);
  while (isdigit((int)*si)) si++;
  skipwhitespace(si);
  V->topseqlen= atol(si);
  /* fprintf(stderr,"Phylip-ileaf: topnseq=%d  topseqlen=%d\n",V->topnseq, V->topseqlen); */

  do {
    getline1(V);
    V->done = feof(V->f);

    if (V->done && !(*V->s)) break;
    si= V->s;
    skipwhitespace(si);
    if (*si != 0) {

      if (first) {  /* collect seq names + seq, as fprintf(outf,"%-10s  ",seqname); */
        (V->nseq)++;
        if (V->nseq >= V->topnseq) first= false;
        sj= V->s+10;  /* past name, start of data */
        if (V->choice == kListSequences) {
          *sj= 0;
          addinfo( si, V);
          }
        else if (V->nseq == V->choice) {
          addseq(sj, V);
          *sj= 0;
          strcpy(V->seqid, si);
          }
        }
      else if ( iline % V->nseq == V->choice -1 ) {
        addseq(si, V);
        }
      iline++;
    }
  } while (!V->done);

  V->allDone = true;
} /*readPhylipInterleaved*/



Local boolean endPhylipSequential( boolean *addend, boolean *ungetend, struct ReadSeqVars *V)
{
  *addend = false;
  *ungetend= false;
  countseq( V->s, V);
  return V->seqlencount >= V->topseqlen;
}

Local void readPhylipSequential(struct ReadSeqVars *V)
{
  short  i;
  char  *si;
  /* sscanf( V->s, "%d%d", &V->topnseq, &V->topseqlen); < ? bad sscan ? */
  si= V->s;
  skipwhitespace(si);
  V->topnseq= atoi(si);
  while (isdigit((int)*si)) si++;
  skipwhitespace(si);
  V->topseqlen= atol(si);
  getline1(V);
  while (!V->allDone) {
    V->seqlencount= 0;
    strncpy(V->seqid, (V->s), 10);
    V->seqid[10]= 0;
    for (i=0; i<10 && V->s[i]; i++) V->s[i]= ' ';
    readLoop(0, true, endPhylipSequential, V);
    if (feof(V->f)) V->allDone = true;
    }
}




Local void readSeqMain(
      struct ReadSeqVars *V,
      const long  skiplines_,
      const short format_)
{
#define tolowerstr(s) { long Itlwr, Ntlwr= strlen(s); \
  for (Itlwr=0; Itlwr<Ntlwr; Itlwr++) s[Itlwr]= to_lower(s[Itlwr]); }

  boolean gotuw;
  long l;

  V->linestart= 0;
  V->matchchar= 0;
  if (V->f == NULL)
    V->err = eFileNotFound;
  else {

    for (l = skiplines_; l > 0; l--) getline1( V);

    do {
      getline1( V);
      for (l= strlen(V->s); (l > 0) && (V->s[l] == ' '); l--) ;
    } while ((l == 0) && !feof(V->f));

    if (feof(V->f)) V->err = eNoData;
    else switch (format_) {
      case kPlain : readPlain(V); break;
      case kIG    : readIG(V); break;
      case kStrider: readStrider(V); break;
      case kGenBank: readGenBank(V); break;
      case kPIR   : readPIR(V); break;
      case kNBRF  : readNBRF(V); break;
      case kPearson: readPearson(V); break;
      case kEMBL  : readEMBL(V); break;
      case kZuker : readZuker(V); break;
      case kOlsen : readOlsen(V); break;
      case kMSF   : readMSF(V); break;

      case kPAUP    : {
        boolean done= false;
        boolean interleaved= false;
        char  *cp;
        /* rewind(V->f); V->nseq= 0; ?? assume it is at top ?? skiplines ... */
        while (!done) {
          getline1( V);
          tolowerstr( V->s);
          if (strstr( V->s, "matrix")) done= true;
          if (strstr( V->s, "interleav")) interleaved= true;
          if (NULL != (cp=strstr( V->s, "ntax=")) )  V->topnseq= atoi(cp+5);
          if (NULL != (cp=strstr( V->s, "nchar=")) )  V->topseqlen= atoi(cp+6);
          if (NULL != (cp=strstr( V->s, "matchchar=")) )  {
            cp += 10;
            if (*cp=='\'') cp++;
            else if (*cp=='"') cp++;
            V->matchchar= *cp;
            }
          }
        if (interleaved) readPAUPinterleaved(V);
        else readPAUPsequential(V);
        }
        break;

      /* kPhylip: ! can't determine in middle of file which type it is...*/
      /* test for interleave or sequential and use Phylip4(ileave) or Phylip2 */
      case kPhylip2:
        readPhylipSequential(V);
        break;
      case kPhylip4: /* == kPhylip3 */
        readPhylipInterleaved(V);
        break;

      default:
        V->err = eUnknownFormat;
        break;

      case kFitch :
        strcpy(V->seqid, V->s); getline1(V);
        readFitch(V);
        break;

      case kGCG:
        do {
          gotuw = (strstr(V->s,"..") != NULL);
          if (gotuw) readUWGCG(V);
          getline1(V);
        } while (!(feof(V->f) || V->allDone));
        break;
      }
    }

  V->filestart= false;
  V->seq[V->seqlen] = 0; /* stick a string terminator on it */
}


char *readSeqFp(
      const short whichEntry_,  /* index to sequence in file */
      FILE  *fp_,   /* pointer to open seq file */
      const long  skiplines_,
      const short format_,      /* sequence file format */
      long  *seqlen_,     /* return seq size */
      short *nseq_,       /* number of seqs in file, for listSeqs() */
      short *error_,      /* return error */
      char  *seqid_)      /* return seq name/info */
{
  struct ReadSeqVars V;

  if (format_ < kMinFormat || format_ > kMaxFormat) {
    *error_ = eUnknownFormat;
    *seqlen_ = 0;
    return NULL;
    }

  V.choice = whichEntry_;
  V.fname  = NULL;  /* don't know */
  V.seq    = (char*) calloc(1, kStartLength+1);
  V.maxseq = kStartLength;
  V.seqlen = 0;
  V.seqid  = seqid_;

  V.f = fp_;
  V.filestart= (ftell( fp_) == 0);
  /* !! in sequential read, must remove current seq position from choice/whichEntry_ counter !! ... */
  if (V.filestart)  V.nseq = 0;
  else V.nseq= *nseq_;  /* track where we are in file...*/

  *V.seqid = '\0';
  V.err = 0;
  V.nseq = 0;
  V.isseqchar = isSeqChar;
  if (V.choice == kListSequences) ; /* leave as is */
  else if (V.choice <= 0) V.choice = 1; /* default ?? */
  V.addit = (V.choice > 0);
  V.allDone = false;

  readSeqMain(&V, skiplines_, format_);

  *error_ = V.err;
  *seqlen_ = V.seqlen;
  *nseq_ = V.nseq;
  return V.seq;
}

char *readSeq(
      const short whichEntry_,  /* index to sequence in file */
      const char  *filename_,   /* file name */
      const long  skiplines_,
      const short format_,      /* sequence file format */
      long  *seqlen_,     /* return seq size */
      short *nseq_,       /* number of seqs in file, for listSeqs() */
      short *error_,      /* return error */
      char  *seqid_)      /* return seq name/info */
{
  struct ReadSeqVars V;

  if (format_ < kMinFormat || format_ > kMaxFormat) {
    *error_ = eUnknownFormat;
    *seqlen_ = 0;
    return NULL;
    }

  V.choice = whichEntry_;
  V.fname  = filename_;  /* don't need to copy string, just ptr to it */
  V.seq    = (char*) calloc(1, kStartLength+1);
  V.maxseq = kStartLength;
  V.seqlen = 0;
  V.seqid  = seqid_;

  V.f = NULL;
  *V.seqid = '\0';
  V.err = 0;
  V.nseq = 0;
  V.isseqchar = isSeqChar;
  if (V.choice == kListSequences) ; /* leave as is */
  else if (V.choice <= 0) V.choice = 1; /* default ?? */
  V.addit = (V.choice > 0);
  V.allDone = false;

  V.f = fopen(V.fname, "r");
  V.filestart= true;

  readSeqMain(&V, skiplines_, format_);

  if (V.f != NULL) fclose(V.f);
  *error_ = V.err;
  *seqlen_ = V.seqlen;
  *nseq_ = V.nseq;
  return V.seq;
}





char *listSeqs(
      const char  *filename_,   /* file name */
      const long skiplines_,
      const short format_,      /* sequence file format */
      short *nseq_,       /* number of seqs in file, for listSeqs() */
      short *error_)      /* return error */
{
  char  seqid[MAXLINE];
  long  seqlen;

  return readSeq( kListSequences, filename_, skiplines_, format_,
                  &seqlen, nseq_, error_, seqid);
}




short seqFileFormat(    /* return sequence format number, see ureadseq.h */
    const char *filename,
    long  *skiplines,   /* return #lines to skip any junk like mail header */
    short *error)       /* return any error value or 0 */
{
  FILE      *fseq;
  short      format;

  fseq  = fopen(filename, "r");
  format= seqFileFormatFp( fseq, skiplines, error);
  if (fseq!=NULL) fclose(fseq);
  return format;
}

short seqFileFormatFp(
    FILE *fseq,
    long  *skiplines,   /* return #lines to skip any junk like mail header */
    short *error)       /* return any error value or 0 */
{
  boolean   foundIG= false, foundStrider= false,
            foundGB= false, foundPIR= false, foundEMBL= false, foundNBRF= false,
            foundPearson= false, foundFitch= false, foundPhylip= false, foundZuker= false,
            gotolsen= false, gotpaup = false, gotasn1 = false, gotuw= false, gotMSF= false,
            isfitch= false,  isphylip= false, done= false;
  short     format= kUnknown;
  int       nlines= 0, k=0, splen= 0, otherlines= 0, aminolines= 0, dnalines= 0;
  char      sp[MAXLINE];
  long      linestart=0;
  int     maxlines2check=5000;

#define ReadOneLine(sp)   \
  { done |= (feof(fseq)); \
    readline( fseq, sp, &linestart);  \
    if (!done) { splen = strlen(sp); ++nlines; } }

  *skiplines = 0;
  *error = 0;
  if (fseq == NULL) { *error = eFileNotFound;  return kNoformat; }

  while ( !done ) {
    ReadOneLine(sp);

    /* check for mailer head & skip past if found */
    if (nlines < 4 && !done) {
      if ((strstr(sp,"From ") == sp) || (strstr(sp,"Received:") == sp)) {
        do {
          /* skip all lines until find one blank line */
          ReadOneLine(sp);
          if (!done) for (k=0; (k<splen) && (sp[k]==' '); k++) ;
          } while ((!done) && (k < splen));
        *skiplines = nlines; /* !? do we want #lines or #bytes ?? */
        }
      }

    if (sp==NULL || *sp==0)
      ; /* nada */

    /* high probability identities: */

    else if ( strstr(sp,"MSF:") && strstr(sp,"Type:") && strstr(sp,"Check:") )
      gotMSF= true;

    else if ((strstr(sp,"..") != NULL) && (strstr(sp,"Check:") != NULL))
      gotuw= true;

    else if (strstr(sp,"identity:   Data:") != NULL)
      gotolsen= true;

    else if ( strstr(sp,"::=") &&
      (strstr(sp,"Bioseq") ||       /* Bioseq or Bioseq-set */
       strstr(sp,"Seq-entry") ||
       strstr(sp,"Seq-submit") ) )  /* can we read submit format? */
          gotasn1= true;

    else if ( strstr(sp,"#NEXUS") == sp )
      gotpaup= true;

    /* uncertain identities: */

    else if (*sp ==';') {
      if (strstr(sp,"Strider") !=NULL) foundStrider= true;
      else foundIG= true;
      }

    else if (strstr(sp,"LOCUS") == sp)
      foundGB= true;
    else if (strstr(sp,"ORIGIN") == sp)
      foundGB= true;

    else if (strstr(sp,"ENTRY   ") == sp)  /* ? also (strcmp(sp,"\\\\\\")==0) */
      foundPIR= true;
    else if (strstr(sp,"SEQUENCE") == sp)
      foundPIR= true;

    else if (*sp == '>') {
      if (sp[3] == ';') foundNBRF= true;
      else foundPearson= true;
      }

    else if (strstr(sp,"ID   ") == sp)
      foundEMBL= true;
    else if (strstr(sp,"SQ   ") == sp)
      foundEMBL= true;

    else if (*sp == '(')
      foundZuker= true;

    else {
      if (nlines - *skiplines == 1) {
        int ispp= 0, ilen= 0;
        sscanf( sp, "%d%d", &ispp, &ilen);
        if (ispp > 0 && ilen > 0) isphylip= true;
        }
      else if (isphylip && nlines - *skiplines == 2) {
        int  tseq;
        tseq= getseqtype(sp+10, strlen(sp+10));
        if ( isalpha((int)*sp)     /* 1st letter in 2nd line must be of a name */
         && (tseq != kOtherSeq))  /* sequence section must be okay */
            foundPhylip= true;
        }

      for (k=0, isfitch= true; isfitch & (k < splen); k++) {
        if (k % 4 == 0) isfitch &= (sp[k] == ' ');
        else isfitch &= (sp[k] != ' ');
        }
      if (isfitch & (splen > 20)) foundFitch= true;

      /* kRNA && kDNA are fairly certain...*/
      switch (getseqtype( sp, splen)) {
        case kOtherSeq: otherlines++; break;
        case kAmino   : if (splen>20) aminolines++; break;
        case kDNA     :
        case kRNA     : if (splen>20) dnalines++; break;
        case kNucleic : break; /* not much info ? */
        }

      }

                    /* pretty certain */
    if (gotolsen) {
      format= kOlsen;
      done= true;
      }
    else if (gotMSF) {
      format= kMSF;
      done= true;
      }
    else if (gotasn1) {
      /* !! we need to look further and return  kASNseqentry | kASNseqset */
      /*
        seqentry key is Seq-entry ::=
        seqset key is Bioseq-set ::=
        ?? can't read these yet w/ ncbi tools ??
          Seq-submit ::=
          Bioseq ::=  << fails both bioseq-seq and seq-entry parsers !
      */
      if (strstr(sp,"Bioseq-set")) format= kASNseqset;
      else if (strstr(sp,"Seq-entry")) format= kASNseqentry;
      else format= kASN1;  /* other form, we can't yet read... */
      done= true;
      }
    else if (gotpaup) {
      format= kPAUP;
      done= true;
      }

    else if (gotuw) {
      if (foundIG) format= kIG;  /* a TOIG file from GCG for certain */
      else format= kGCG;
      done= true;
      }

    else if ((dnalines > 1) || done || (nlines > maxlines2check)) {
          /* decide on most likely format */
          /* multichar idents: */
      if (foundStrider) format= kStrider;
      else if (foundGB) format= kGenBank;
      else if (foundPIR) format= kPIR;
      else if (foundEMBL) format= kEMBL;
      else if (foundNBRF) format= kNBRF;
          /* single char idents: */
      else if (foundIG) format= kIG;
      else if (foundPearson) format= kPearson;
      else if (foundZuker) format= kZuker;
          /* digit ident: */
      else if (foundPhylip) format= kPhylip;
          /* spacing ident: */
      else if (foundFitch) format= kFitch;
          /* no format chars: */
      else if (otherlines > 0) format= kUnknown;
      else if (dnalines > 1) format= kPlain;
      else if (aminolines > 1) format= kPlain;
      else format= kUnknown;

      done= true;
      }

    /* need this for possible long header in olsen format */
     else if (strstr(sp,"): ") != NULL)
       maxlines2check++;
    }

  if (format == kPhylip) {
    /* check for interleaved or sequential -- really messy */
    int tname, tseq;
    long i, j, nspp= 0, nlen= 0, ilen, leaf= 0, seq= 0;
    char  *ps;

    rewind(fseq);
    for (i=0; i < *skiplines; i++) ReadOneLine(sp);
    nlines= 0;
    ReadOneLine(sp);
    sscanf( sp, "%ld%ld", &nspp, &nlen);
    ReadOneLine(sp); /* 1st seq line */
    for (ps= sp+10, ilen=0; *ps!=0; ps++) if (isprint((int)*ps)) ilen++;

    for (i= 1; i<nspp; i++) {
      ReadOneLine(sp);

      tseq= getseqtype(sp+10, strlen(sp+10));
      tname= getseqtype(sp, 10);
      for (j=0, ps= sp; isspace((int)*ps) && j<10; ps++, j++);
      for (ps= sp; *ps!=0; ps++) if (isprint((int)*ps)) ilen++;

      /* find probable interleaf or sequential ... */
      if (j>=9) seq += 10; /* pretty certain not ileaf */
      else {
        if (tseq != tname) leaf++; else seq++;
        if (tname == kDNA || tname == kRNA) seq++; else leaf++;
        }

      if (ilen <= nlen && j<9) {
        if (tname == kOtherSeq) leaf += 10;
        else if (tname == kAmino || tname == kDNA || tname == kRNA) seq++; else leaf++;
        }
      else if (ilen > nlen) {
        ilen= 0;
        }
      }
    for ( nspp *= 2 ; i<nspp; i++) {  /* this should be only bases if interleaf */
      ReadOneLine(sp);

      tseq= getseqtype(sp+10, strlen(sp+10));
      tname= getseqtype(sp, 10);
      for (ps= sp; *ps!=0; ps++) if (isprint((int)*ps)) ilen++;
      for (j=0, ps= sp; isspace((int)*ps) && j<10; ps++, j++);
      if (j<9) {
        if (tname == kOtherSeq) seq += 10;
        if (tseq != tname) seq++; else leaf++;
        if (tname == kDNA || tname == kRNA) leaf++; else seq++;
        }
      if (ilen > nlen) {
        if (j>9) leaf += 10; /* must be a name here for sequent */
        else if (tname == kOtherSeq) seq += 10;
        ilen= 0;
        }
      }

    if (leaf > seq) format= kPhylip4;
    else format= kPhylip2;
    }

  return(format);
#undef  ReadOneLine
} /* SeqFileFormat */




unsigned long GCGchecksum( const char *seq, const long seqlen, unsigned long *checktotal)
/* GCGchecksum */
{
  register long  i, check = 0, count = 0;

  for (i = 0; i < seqlen; i++) {
    count++;
    check += count * to_upper(seq[i]);
    if (count == 57) count = 0;
    }
  check %= 10000;
  *checktotal += check;
  *checktotal %= 10000;
  return check;
}

/* Table of CRC-32's of all single byte values (made by makecrc.c of ZIP source) */
const unsigned long crctab[] = {
  0x00000000UL, 0x77073096UL, 0xee0e612cUL, 0x990951baUL, 0x076dc419UL,
  0x706af48fUL, 0xe963a535UL, 0x9e6495a3UL, 0x0edb8832UL, 0x79dcb8a4UL,
  0xe0d5e91eUL, 0x97d2d988UL, 0x09b64c2bUL, 0x7eb17cbdUL, 0xe7b82d07UL,
  0x90bf1d91UL, 0x1db71064UL, 0x6ab020f2UL, 0xf3b97148UL, 0x84be41deUL,
  0x1adad47dUL, 0x6ddde4ebUL, 0xf4d4b551UL, 0x83d385c7UL, 0x136c9856UL,
  0x646ba8c0UL, 0xfd62f97aUL, 0x8a65c9ecUL, 0x14015c4fUL, 0x63066cd9UL,
  0xfa0f3d63UL, 0x8d080df5UL, 0x3b6e20c8UL, 0x4c69105eUL, 0xd56041e4UL,
  0xa2677172UL, 0x3c03e4d1UL, 0x4b04d447UL, 0xd20d85fdUL, 0xa50ab56bUL,
  0x35b5a8faUL, 0x42b2986cUL, 0xdbbbc9d6UL, 0xacbcf940UL, 0x32d86ce3UL,
  0x45df5c75UL, 0xdcd60dcfUL, 0xabd13d59UL, 0x26d930acUL, 0x51de003aUL,
  0xc8d75180UL, 0xbfd06116UL, 0x21b4f4b5UL, 0x56b3c423UL, 0xcfba9599UL,
  0xb8bda50fUL, 0x2802b89eUL, 0x5f058808UL, 0xc60cd9b2UL, 0xb10be924UL,
  0x2f6f7c87UL, 0x58684c11UL, 0xc1611dabUL, 0xb6662d3dUL, 0x76dc4190UL,
  0x01db7106UL, 0x98d220bcUL, 0xefd5102aUL, 0x71b18589UL, 0x06b6b51fUL,
  0x9fbfe4a5UL, 0xe8b8d433UL, 0x7807c9a2UL, 0x0f00f934UL, 0x9609a88eUL,
  0xe10e9818UL, 0x7f6a0dbbUL, 0x086d3d2dUL, 0x91646c97UL, 0xe6635c01UL,
  0x6b6b51f4UL, 0x1c6c6162UL, 0x856530d8UL, 0xf262004eUL, 0x6c0695edUL,
  0x1b01a57bUL, 0x8208f4c1UL, 0xf50fc457UL, 0x65b0d9c6UL, 0x12b7e950UL,
  0x8bbeb8eaUL, 0xfcb9887cUL, 0x62dd1ddfUL, 0x15da2d49UL, 0x8cd37cf3UL,
  0xfbd44c65UL, 0x4db26158UL, 0x3ab551ceUL, 0xa3bc0074UL, 0xd4bb30e2UL,
  0x4adfa541UL, 0x3dd895d7UL, 0xa4d1c46dUL, 0xd3d6f4fbUL, 0x4369e96aUL,
  0x346ed9fcUL, 0xad678846UL, 0xda60b8d0UL, 0x44042d73UL, 0x33031de5UL,
  0xaa0a4c5fUL, 0xdd0d7cc9UL, 0x5005713cUL, 0x270241aaUL, 0xbe0b1010UL,
  0xc90c2086UL, 0x5768b525UL, 0x206f85b3UL, 0xb966d409UL, 0xce61e49fUL,
  0x5edef90eUL, 0x29d9c998UL, 0xb0d09822UL, 0xc7d7a8b4UL, 0x59b33d17UL,
  0x2eb40d81UL, 0xb7bd5c3bUL, 0xc0ba6cadUL, 0xedb88320UL, 0x9abfb3b6UL,
  0x03b6e20cUL, 0x74b1d29aUL, 0xead54739UL, 0x9dd277afUL, 0x04db2615UL,
  0x73dc1683UL, 0xe3630b12UL, 0x94643b84UL, 0x0d6d6a3eUL, 0x7a6a5aa8UL,
  0xe40ecf0bUL, 0x9309ff9dUL, 0x0a00ae27UL, 0x7d079eb1UL, 0xf00f9344UL,
  0x8708a3d2UL, 0x1e01f268UL, 0x6906c2feUL, 0xf762575dUL, 0x806567cbUL,
  0x196c3671UL, 0x6e6b06e7UL, 0xfed41b76UL, 0x89d32be0UL, 0x10da7a5aUL,
  0x67dd4accUL, 0xf9b9df6fUL, 0x8ebeeff9UL, 0x17b7be43UL, 0x60b08ed5UL,
  0xd6d6a3e8UL, 0xa1d1937eUL, 0x38d8c2c4UL, 0x4fdff252UL, 0xd1bb67f1UL,
  0xa6bc5767UL, 0x3fb506ddUL, 0x48b2364bUL, 0xd80d2bdaUL, 0xaf0a1b4cUL,
  0x36034af6UL, 0x41047a60UL, 0xdf60efc3UL, 0xa867df55UL, 0x316e8eefUL,
  0x4669be79UL, 0xcb61b38cUL, 0xbc66831aUL, 0x256fd2a0UL, 0x5268e236UL,
  0xcc0c7795UL, 0xbb0b4703UL, 0x220216b9UL, 0x5505262fUL, 0xc5ba3bbeUL,
  0xb2bd0b28UL, 0x2bb45a92UL, 0x5cb36a04UL, 0xc2d7ffa7UL, 0xb5d0cf31UL,
  0x2cd99e8bUL, 0x5bdeae1dUL, 0x9b64c2b0UL, 0xec63f226UL, 0x756aa39cUL,
  0x026d930aUL, 0x9c0906a9UL, 0xeb0e363fUL, 0x72076785UL, 0x05005713UL,
  0x95bf4a82UL, 0xe2b87a14UL, 0x7bb12baeUL, 0x0cb61b38UL, 0x92d28e9bUL,
  0xe5d5be0dUL, 0x7cdcefb7UL, 0x0bdbdf21UL, 0x86d3d2d4UL, 0xf1d4e242UL,
  0x68ddb3f8UL, 0x1fda836eUL, 0x81be16cdUL, 0xf6b9265bUL, 0x6fb077e1UL,
  0x18b74777UL, 0x88085ae6UL, 0xff0f6a70UL, 0x66063bcaUL, 0x11010b5cUL,
  0x8f659effUL, 0xf862ae69UL, 0x616bffd3UL, 0x166ccf45UL, 0xa00ae278UL,
  0xd70dd2eeUL, 0x4e048354UL, 0x3903b3c2UL, 0xa7672661UL, 0xd06016f7UL,
  0x4969474dUL, 0x3e6e77dbUL, 0xaed16a4aUL, 0xd9d65adcUL, 0x40df0b66UL,
  0x37d83bf0UL, 0xa9bcae53UL, 0xdebb9ec5UL, 0x47b2cf7fUL, 0x30b5ffe9UL,
  0xbdbdf21cUL, 0xcabac28aUL, 0x53b39330UL, 0x24b4a3a6UL, 0xbad03605UL,
  0xcdd70693UL, 0x54de5729UL, 0x23d967bfUL, 0xb3667a2eUL, 0xc4614ab8UL,
  0x5d681b02UL, 0x2a6f2b94UL, 0xb40bbe37UL, 0xc30c8ea1UL, 0x5a05df1bUL,
  0x2d02ef8dL
};

unsigned long CRC32checksum(const char *seq, const long seqlen, unsigned long *checktotal)
/*CRC32checksum: modified from CRC-32 algorithm found in ZIP compression source */
{
  register unsigned long c = 0xffffffffUL;
  register long n = seqlen;

  /* tbailey; 11-11-96: move seq++ out of macro! */
  while (n--) {
    c = crctab[((int)c ^ (to_upper(*seq))) & 0xff] ^ (c >> 8); 
    seq++;
  }
  c= c ^ 0xffffffffUL;
  *checktotal += c;
  return c;
}




short getseqtype( const char *seq, const long seqlen)
{ /* return sequence kind: kDNA, kRNA, kProtein, kOtherSeq, ??? */
  char  c;
  short i, maxtest;
  short na = 0, aa = 0, po = 0, nt = 0, nu = 0, ns = 0, no = 0;

  maxtest = min(300, seqlen);
  for (i = 0; i < maxtest; i++) {
    c = to_upper(seq[i]);
    if (strchr(protonly, c)) po++;
    else if (strchr(primenuc,c)) {
      na++;
      if (c == 'T') nt++;
      else if (c == 'U') nu++;
      }
    else if (strchr(aminos,c)) aa++;
    else if (strchr(seqsymbols,c)) ns++;
    else if (isalpha((int)c)) no++;
    }

  if ((no > 0) || (po+aa+na == 0)) return kOtherSeq;
  /* ?? test for probability of kOtherSeq ?, e.g.,
  else if (po+aa+na / maxtest < 0.70) return kOtherSeq;
  */
  else if (po > 0) return kAmino;
  else if (aa == 0) {
    if (nu > nt) return kRNA;
    else return kDNA;
    }
  else if (na > aa) return kNucleic;
  else return kAmino;
} /* getseqtype */


char* compressSeq( const char gapc, const char *seq, const long seqlen, long *newlen)
{
  register char *a, *b;
  register long i;
  char  *newseq;

  *newlen= 0;
  if (!seq) return NULL;
  newseq = (char*) malloc(seqlen+1);
  if (!newseq) return NULL;
  for (a= (char*)seq, b=newseq, i=0; *a!=0; a++)
    if (*a != gapc) {
      *b++= *a;
      i++;
      }
  *b= '\0';
  newseq = (char*) realloc(newseq, i+1);
  *newlen= i;
  return newseq;
}



/***
char *rtfhead = "{\\rtf1\\defformat\\mac\\deff2 \
{\\fonttbl\
  {\\f1\\fmodern Courier;}{\\f2\\fmodern Monaco;}\
  {\\f3\\fswiss Helvetica;}{\\f4\\fswiss Geneva;}\
  {\\f5\\froman Times;}{\\f6\\froman Palatino;}\
  {\\f7\\froman New Century Schlbk;}{\\f8\\ftech Symbol;}}\
{\\stylesheet\
  {\\s1 \\f5\\fs20 \\sbasedon0\\snext1 name;}\
  {\\s2 \\f3\\fs20 \\sbasedon0\\snext2 num;}\
  {\\s3 \\f1\\f21 \\sbasedon0\\snext3 seq;}}";

char *rtftail = "}";
****/

short writeSeq(FILE *outf, const char *seq, const long seqlen,
                const short outform, const char *seqid)
/* dump sequence to standard output */
{
  const short kSpaceAll = -9;
#define kMaxseqwidth  250

  boolean baseonlynum= false; /* nocountsymbols -- only count true bases, not "-" */
  short  numline = 0; /* only true if we are writing seq number line (for interleave) */
  boolean numright = false, numleft = false;
  boolean nameright = false, nameleft = false;
  short   namewidth = 8, numwidth = 8;
  short   spacer = 0, width  = 50, tab = 0;
  /* new parameters: width, spacer, those above... */

  short linesout = 0, seqtype = kNucleic;
  long  i, j, l, l1, ibase;
  char  idword[31], endstr[14];
  char  seqnamestore[MAXLINE], *seqname = seqnamestore;
  char  s[kMaxseqwidth], *cp=NULL;
  char  nameform[10], numform[10], nocountsymbols[10];
  unsigned long checksum = 0, checktotal = 0;

  gPretty.atseq++;
  skipwhitespace(seqid);
  /*l = min(128, strlen(seqid)); */
  strlcpy(seqnamestore, seqid, MAXLINE);
  /*seqname[l] = 0; */

  sscanf( seqname, "%30s", idword);
  sprintf(numform, "%ld", seqlen);
  numwidth= strlen(numform)+1;
  nameform[0]= '\0';

  if (strstr(seqname,"checksum") != NULL) {
    cp = strstr(seqname,"bases");
    if (cp!=NULL) {
      for ( ; (cp!=seqname) && (*cp!=','); cp--) ;
      if (cp!=seqname) *cp=0;
      }
    }

  strcpy( endstr,"");
  l1 = 0;

  if (outform == kGCG || outform == kMSF)
    checksum = GCGchecksum(seq, seqlen, &checktotal);
  else
    checksum = seqchecksum(seq, seqlen, &checktotal);

  switch (outform) {

    case kPlain:
    case kUnknown:    /* no header, just sequence */
      strcpy(endstr,"\n"); /* end w/ extra blank line */
      break;

    case kOlsen:  /* Olsen seq. editor takes plain nucs OR Genbank  */
    case kGenBank:
      fprintf(outf,"LOCUS       %s       %ld bp\n", idword, seqlen);
      fprintf(outf,"DEFINITION  %s, %ld bases, %X checksum.\n", 
        seqname, seqlen, (int)checksum);
   /* fprintf(outf,"ACCESSION   %s\n", accnum); */
      fprintf(outf,"ORIGIN      \n");
      spacer = 11;
      numleft = true;
      numwidth = 8;  /* dgg. 1Feb93, patch for GDE fail to read short numwidth */
      strcpy(endstr, "\n//");
      linesout += 4;
      break;

    case kPIR:
      /* somewhat like genbank... \\\*/
      /* fprintf(outf,"\\\\\\\n"); << only at top of file, not each entry... */
      fprintf(outf,"ENTRY           %s \n", idword);
      fprintf(outf,"TITLE           %s, %ld bases, %X checksum.\n", 
        seqname, seqlen, (int)checksum);
   /* fprintf(outf,"ACCESSION       %s\n", accnum); */
      fprintf(outf,"SEQUENCE        \n");
      numwidth = 7;
      width= 30;
      spacer = kSpaceAll;
      numleft = true;
      strcpy(endstr, "\n///");
      /* run a top number line for PIR */
      for (j=0; j<numwidth; j++) fputc(' ',outf);
      for (j= 5; j<=width; j += 5) fprintf(outf,"%10ld",j);
      fputc('\n',outf);
      linesout += 5;
      break;

    case kNBRF:
      if (getseqtype(seq, seqlen) == kAmino)
        fprintf(outf,">P1;%s\n", idword);
      else
        fprintf(outf,">DL;%s\n", idword);
      fprintf(outf,"%s, %ld bases, %X checksum.\n", 
        seqname, seqlen, (int)checksum);
      spacer = 11;
      strcpy(endstr,"*\n");
      linesout += 3;
      break;

    case kEMBL:
      fprintf(outf,"ID   %s\n", idword);
  /*  fprintf(outf,"AC   %s\n", accnum); */
      fprintf(outf,"DE   %s, %ld bases, %X checksum.\n", 
        seqname, seqlen, (int)checksum);
      fprintf(outf,"SQ             %ld BP\n", seqlen);
      strcpy(endstr, "\n//"); /* 11Oct90: bug fix*/
      tab = 4;     /** added 31jan91 */
      spacer = 11; /** added 31jan91 */
      width = 60;
      linesout += 4;
      break;

    case kGCG:
      fprintf(outf,"%s\n", seqname);
   /* fprintf(outf,"ACCESSION   %s\n", accnum); */
      fprintf(outf,"    %s  Length: %ld  (today)  Check: %d  ..\n", 
        idword, seqlen, (int)checksum);
      spacer = 11;
      numleft = true;
      strcpy(endstr, "\n");  /* this is insurance to help prevent misreads at eof */
      linesout += 3;
      break;

    case kStrider: /* ?? map ?*/
      fprintf(outf,"; ### from DNA Strider ;-)\n");
      fprintf(outf,"; DNA sequence  %s, %ld bases, %X checksum.\n;\n", 
        seqname, seqlen, (int)checksum);
      strcpy(endstr, "\n//");
      linesout += 3;
      break;

    case kFitch:
      fprintf(outf,"%s, %ld bases, %X checksum.\n", 
        seqname, seqlen, (int)checksum);
      spacer = 4;
      width = 60;
      linesout += 1;
      break;

    case kPhylip2:
    case kPhylip4:
      /* this is version 3.2/3.4 -- simplest way to write
        version 3.3 is to write as version 3.2, then
        re-read file and interleave the species lines */
      if (strlen(idword)>10) idword[10] = 0;
      fprintf(outf,"%-10s  ",idword);
      l1  = -1;
      tab = 12;
      spacer = 11;
      break;

    case kASN1:
      seqtype= getseqtype(seq, seqlen);
      switch (seqtype) {
        case kDNA     : cp= "dna"; break;
        case kRNA     : cp= "rna"; break;
        case kNucleic : cp= "na"; break;
        case kAmino   : cp= "aa"; break;
        case kOtherSeq: cp= "not-set"; break;
        }
      fprintf(outf,"  seq {\n");
      fprintf(outf,"    id { local id %d },\n", gPretty.atseq);
      fprintf(outf,"    descr { title \"%s\" },\n", seqid);
      fprintf(outf,"    inst {\n");
      fprintf(outf,"      repr raw, mol %s, length %ld,\ntopology linear,\n", cp, seqlen);
      fprintf(outf,"      seq-data\n");
      if (seqtype == kAmino)
        fprintf(outf,"        iupacaa \"");
      else
        fprintf(outf,"        iupacna \"");
      l1  = 17;
      spacer = 0;
      width  = 78;
      tab  = 0;
      strcpy(endstr,"\"\n      } } ,");
      linesout += 7;
      break;

    case kPAUP:
      nameleft= true;
      namewidth = 9;
      spacer = 21;
      width  = 100;
      tab  = 0; /* 1; */
      /* strcpy(endstr,";\nend;"); << this is end of all seqs.. */
      /* do a header comment line for paup */
      fprintf(outf,"[Name: %-16s  Len:%6ld  Check: %8X]\n", 
        idword, seqlen, (int)checksum);
      linesout += 1;
      break;

    case kPretty:
      numline= gPretty.numline;
      baseonlynum= gPretty.baseonlynum;
      namewidth = gPretty.namewidth;
      numright = gPretty.numright;
      numleft = gPretty.numleft;
      nameright = gPretty.nameright;
      nameleft = gPretty.nameleft;
      spacer = gPretty.spacer + 1;
      width  = gPretty.seqwidth;
      tab  = gPretty.tab;
      /* also add rtf formatting w/ font, size, style */
      if (gPretty.nametop) {
        fprintf(outf,"Name: %-16s  Len:%6ld  Check: %8X\n", 
         idword, seqlen, (int)checksum);
        linesout++;
        }
      break;

    case kMSF:
      fprintf(outf," Name: %-16s Len:%6ld  Check: %5d  Weight:  1.00\n",
                    idword, seqlen, (int)checksum);
      linesout++;
      nameleft= true;
      namewidth= 15; /* need MAX namewidth here... */
      sprintf(nameform, "%%+%ds ",namewidth);
      spacer = 11;
      width  = 50;
      tab = 0; /* 1; */
      break;

    case kIG:
      fprintf(outf,";%s, %ld bases, %X checksum.\n", 
        seqname, seqlen, (int)checksum);
      fprintf(outf,"%s\n", idword);
      strcpy(endstr,"1"); /* == linear dna */
      linesout += 2;
      break;

    default :
    case kZuker: /* don't attempt Zuker's ftn format */
    case kPearson:
      fprintf(outf,">%s %ld bases, %X checksum.\n", 
        seqname, seqlen, (int)checksum);
      linesout += 1;
      break;
    }

  if (*nameform==0) sprintf(nameform, "%%%d.%ds ",namewidth,namewidth);
  if (numline) sprintf(numform, "%%%ds ",numwidth);
  else sprintf(numform, "%%%dd ",numwidth);
  strcpy( nocountsymbols, kNocountsymbols);
  if (baseonlynum) {
    if (strchr(nocountsymbols,gPretty.gapchar)==NULL) {
      strcat(nocountsymbols," ");
      nocountsymbols[strlen(nocountsymbols)-1]= gPretty.gapchar;
      }
    if (gPretty.domatch && (cp=strchr(nocountsymbols,gPretty.matchchar))!=NULL) {
      *cp= ' ';
      }
    }

  if (numline) {
   *idword= 0;
   }

  width = min(width,kMaxseqwidth);
  for (i=0, l=0, ibase = 1; i < seqlen; ) {

    if (l1 < 0) l1 = 0;
    else if (l1 == 0) {
      if (nameleft) fprintf(outf, nameform, idword);
      if (numleft) { if (numline) fprintf(outf, numform, "");
                    else fprintf(outf, numform, ibase);}
      for (j=0; j<tab; j++) fputc(' ',outf);
      }

    l1++;                 /* don't count spaces for width*/
    if (numline) {
      if (spacer==kSpaceAll || (spacer != 0 && (l+1) % spacer == 1)) {
        if (numline==1) fputc(' ',outf);
        s[l++] = ' ';
        }
      if (l1 % 10 == 1 || l1 == width) {
        if (numline==1) fprintf(outf,"%-9ld ",i+1);
        s[l++]= '|'; /* == put a number here */
        }
      else s[l++]= ' ';
      i++;
      }

    else {
      if (spacer==kSpaceAll || (spacer != 0 && (l+1) % spacer == 1))
        s[l++] = ' ';
      if (!baseonlynum) ibase++;
      else if (0==strchr(nocountsymbols,seq[i])) ibase++;
      s[l++] = seq[i++];
      }

    if (l1 == width || i == seqlen) {
      if (outform==kPretty) for ( ; l1<width; l1++) {
        if (spacer==kSpaceAll || (spacer != 0 && (l+1) % spacer == 1))
          s[l++] = ' ';
        s[l++]=' '; /* pad w/ blanks */
        }
      s[l] = '\0';
      l = 0; l1 = 0;

      if (numline) {
        if (numline==2) fprintf(outf,"%s",s); /* finish numberline ! and | */
        }
      else {
        if (i == seqlen) fprintf(outf,"%s%s",s,endstr);
        else fprintf(outf,"%s",s);
        if (numright || nameright) fputc(' ',outf);
        if (numright)  fprintf(outf,numform, ibase-1);
        if (nameright) fprintf(outf, nameform,idword);
        }
      fputc('\n',outf);
      linesout++;
      }
    }
  return linesout;
}  /*writeSeq*/



/* End file: ureadseq.c */
