/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
   ds.c
   deep-shallow sorting algorithm. these routines are taken mainly by 
   Seward's d_copyEQ_u12.c 
   Buckets are sorted one at a time; when  bucket Y is completely sorted
   (except for YY) "pointer copying" is used to sort small buckets of 
   the form XY, X=A..Z
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */ 
#include <stdio.h>
#include <stdlib.h>
#include "common.h"

// ----- external variables ----------
extern int Anchor_dist;
extern int Shallow_limit;
extern int _ds_Verbose; 

// ----- macros and costants --------------
#define BIGFREQ(b) (ftab[((b)+1) << 8] - ftab[(b) << 8])


/* ------ "local" global variables ------- */
Int32  Text_size;           // size of input string 
UChar  *Text;               // input string+ overshoot
Int32  *Sa;                 // suffix array
UChar  *Upper_text_limit;   // Text+Text_size
Int32  *Anchor_rank;        // rank (in the sorted suffixes of the  
                            // anchor points (-1 if rank is unknown))
UInt16  *Anchor_offset;     // offset (wrt to the anchor) of the suffix
                            // whose rank is in Anchor_rank. 
Int32 Anchor_num;           // number of anchor points
Int32 ftab [65537];      
Int32 runningOrder[256];


/* ------------------------------------------------------------------------
   The use of Anchor_rank[] and Anchor_offset is the following: 
   Anchor_rank[i] is either -1 or contains the rank (in the list of
   the sorted suffixes of the suffix starting at position 
     i*Anchor_dist + Anchor_offset[i].
   Initially Anchor_rank[i] = -1 and Anchor_offset[i]=Anchor_dist,
   then, if a suffix in position t (i*Anchor_dist <= t < (i+1)*Anchor_dist)
   appears to be in a large group which is sorted, the rank of
   t is stored in Anchor_rank[i], and the value t-(i*Anchor_dist) 
   is written to Anchor_offset[i]. Both vaulues can be later updated,
   but the value in Anchor_offset[i] can only decrease, so no further
   changes are done when Anchor_offset[i] is = 0. The invariant is:
   if Anchor_rank[i]>=0 then 
       Sa[Anchor_rank[i]]=i*Anchor_dist+Anchor_offset[i]
   -------------------------------------------------------------------------*/
 

  
// ---------- declare and init counters --------------
Int32 Calls_helped_sort=0;     
Int32 Calls_anchor_sort_forw=0;     
Int32 Calls_anchor_sort_backw=0;    
Int32 Calls_pseudo_anchor_sort_forw=0;      
Int32 Calls_deep_sort=0;     

// --- static prototypes (gcc4 requires they are not inside functions) 
static void check_ordering(int, int);
static void calc_running_order ( void );


/* ************************************************************
   This is the main deep/shallow suffix sorting routines
   It divides the suffixes in buckets according to the 
   first two characters. Some of the buckets are then sorted 
   calling shallow_sort(). When all buckets of kind ax (x!=a) are 
   sorted, we use this ordering to sort all suffixes in the 
   buckets ya (for any y including y=a).
   ************************************************************* */
void ds_ssort(UChar *x, Int32 *p, Int32 n, int maxcmp)
{
  void shallow_sort(Int32 *, int, int);
  int compute_overshoot(void), overshoot;
  Int32  i, j, ss, sb, k;
  UChar  c1, c2;
  Bool   bigDone[256];
  Int32  copyStart[256];
  Int32  copyEnd  [256];
  Int32  numQSorted = 0;

  // ------ set some global variable ------
  Text=x;
  Text_size=n;
  Sa = p;
  Upper_text_limit = Text + Text_size;
  // ------ fill overshoot area
  overshoot = compute_overshoot();
  for(i=n;i<n+overshoot;i++) Text[i]=0; 
  
  // ------ init array containing positions of anchors
  if(Anchor_dist==0) {
    Anchor_num=0; Anchor_rank=NULL; Anchor_offset=NULL;
  }
  else {
    Anchor_num = 2 + (n-1)/Anchor_dist;  // see comment for helped_sort() 
    Anchor_rank = (Int32 *) malloc(Anchor_num*sizeof(Int32));
    Anchor_offset = (UInt16 *) malloc(Anchor_num*sizeof(UInt16));
    if(!Anchor_rank || !Anchor_offset) {
      fprintf(stderr, "malloc failed (ds_sort)\n");
      exit(1);
    }
    for(i=0;i<Anchor_num;i++) {
      Anchor_rank[i]= -1;               // pos of anchors is initially unknown
      Anchor_offset[i] = Anchor_dist;   // maximum possible value
    }
  }

  // ---------- init ftab ------------------
  for (i = 0; i <= 65536; i++) ftab[i] = 0;
  c1 = Text[0];
  for (i = 1; i <= Text_size; i++) {
    c2 = Text[i];
    ftab[(c1 << 8) + c2]++;
    c1 = c2;
  }
  for (i = 1; i <= 65536; i++) ftab[i] += ftab[i-1];

  // -------- sort suffixes considering only the first two chars 
  c1 = Text[0];
  for (i = 0; i < Text_size; i++) {
    c2 = Text[i+1];
    j = (c1 << 8) + c2;
    c1 = c2;
    ftab[j]--;
    Sa[ftab[j]] = i;
  }

  /* decide on the running order */
  calc_running_order();
  for (i = 0; i < 256; i++) bigDone[i] = False;

   /* Really do the suffix sorting */
  for (i = 0; i <= 255; i++) {

    /*--
      Process big buckets, starting with the least full.
      --*/
    ss = runningOrder[i];
    if(_ds_Verbose>2)
      fprintf(stderr,"group %3d;  size %d\n",ss,BIGFREQ(ss)&CLEARMASK); 

    /*--
      Complete the big bucket [ss] by sorting
      any unsorted small buckets [ss, j].  Hopefully
      previous pointer-scanning phases have already
      completed many of the small buckets [ss, j], so
      we don't have to sort them at all.
      --*/
    for (j = 0; j <= 255; j++) {
      if (j != ss) {
	sb = (ss << 8) + j;
	if ( ! (ftab[sb] & SETMASK) ) {
	  Int32 lo = ftab[sb]   & CLEARMASK;
	  Int32 hi = (ftab[sb+1] & CLEARMASK) - 1;
	  if (hi > lo) {
	    if (_ds_Verbose>2)
	      fprintf(stderr,"sorting [%02x, %02x], done %d "
			"this %d\n", ss, j, numQSorted, hi - lo + 1 );
	    shallow_sort(Sa+lo, hi-lo+1,maxcmp);
            #if 0
	    check_ordering(lo, hi);
            #endif
	    numQSorted += ( hi - lo + 1 );
	  }
	}
	ftab[sb] |= SETMASK;
      }
    }
    assert (!bigDone[ss]);
    // ------ now order small buckets of type [xx,ss]  --------
    {
      for (j = 0; j <= 255; j++) {
	copyStart[j] =  ftab[(j << 8) + ss]     & CLEARMASK;
	copyEnd  [j] = (ftab[(j << 8) + ss + 1] & CLEARMASK) - 1;
      }
      // take care of the virtual -1 char in position Text_size+1
      if(ss==0) {
	k=Text_size-1;
	c1 = Text[k];
	if (!bigDone[c1])
	  Sa[ copyStart[c1]++ ] = k;
      }
      for (j = ftab[ss << 8] & CLEARMASK; j < copyStart[ss]; j++) {
	k = Sa[j]-1; if (k < 0) continue;  
	c1 = Text[k];
	if (!bigDone[c1])
	  Sa[ copyStart[c1]++ ] = k;
      }
      for (j = (ftab[(ss+1) << 8] & CLEARMASK) - 1; j > copyEnd[ss]; j--) {
	k = Sa[j]-1; if (k < 0) continue;
	c1 = Text[k];
	if (!bigDone[c1]) 
	  Sa[ copyEnd[c1]-- ] = k;
      }
    }
    assert (copyStart[ss] - 1 == copyEnd[ss]);
    for (j = 0; j <= 255; j++) ftab[(j << 8) + ss] |= SETMASK;
    bigDone[ss] = True;
  }
  if (_ds_Verbose) {
    fprintf(stderr, "\t %d pointers, %d sorted, %d scanned\n",
	      Text_size, numQSorted, Text_size - numQSorted );
    fprintf(stderr, "\t %d calls to helped_sort\n",Calls_helped_sort);      
    fprintf(stderr, "\t %d calls to anchor_sort (forward)\n",
	    Calls_anchor_sort_forw);      
    fprintf(stderr, "\t %d calls to anchor_sort (backward)\n",
	    Calls_anchor_sort_backw);      
    fprintf(stderr, "\t %d calls to pseudo_anchor_sort (forward)\n",
    	    Calls_pseudo_anchor_sort_forw);      
    fprintf(stderr, "\t %d calls to deep_sort\n",Calls_deep_sort);      
  }
  // ---- done! ---------------------------------------- 
  free(Anchor_offset);
  free(Anchor_rank);
}



/* ****************************************************************
   compute running =(sorting) order for big buckets: start with 
   the least full and proceed to the largest one.
   The sorting is done using shellsort
   **************************************************************** */ 
static
void calc_running_order ( void )
{
   Int32 i, j;
   for (i = 0; i <= 255; i++) runningOrder[i] = i;

   {
      Int32 vv;
      Int32 h = 1;
      do h = 3 * h + 1; while (h <= 256);
      do {
         h = h / 3;
         for (i = h; i <= 255; i++) {
            vv = runningOrder[i];
            j = i;
            while ( BIGFREQ(runningOrder[j-h]) > BIGFREQ(vv) ) {
               runningOrder[j] = runningOrder[j-h];
               j = j - h;
               if (j <= (h - 1)) goto zero;
            }
            zero:
            runningOrder[j] = vv;
         }
      } while (h != 1);
   }
}


#if 0
/* *******************************************************
   check that the suffixes in position [lo,hi] are in they
   correct order
   ******************************************************** */
static
void check_ordering(int lo, int hi)
{
  int j1,jj,error;

  error=0;
  for(j1=lo;j1<hi;j1++) {
    if (scmp3(Text+Sa[j1], Text+Sa[j1+1], &jj, 
	      MIN(Text_size-Sa[j1],Text_size-Sa[j1+1]))>=0) {
      for(jj=0;jj<10;jj++) 
	pretty_putchar(Text[Sa[j1]+jj]);
      printf("\n");
      for(jj=0;jj<10;jj++) 
	pretty_putchar(Text[Sa[j1+1]+jj]);
      printf("\n");
      error++;
    }
  }
  if(error>0) {
    printf("----------- start ----------\n");
    for(j1=lo;j1<=hi;j1++) {
      for(jj=0;jj<10;jj++) 
	pretty_putchar(Text[Sa[j1]+jj]);
      printf("\n");
    }
    printf("----------- end ------------\n\n");
  }  
}
#endif











