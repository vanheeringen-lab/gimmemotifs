/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   shallow.c  
   This is the multikey quicksort from bentley-sedgewick modified 
   so that it stops recursion when depth reaches  Shallow_limit 
   (that is when two or more suffixes have Shallow_limit chars in common).
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "common.h"

// ----- external variables ----------
extern UChar *Text;                // start of input string
extern UChar *Upper_text_limit;    // Text+Text_size
extern Int32 _ds_Word_size;        // # of bytes in a word in mkq
extern Int32 Mk_qs_thresh;         // recursion limit for mk quicksort: 
                                   // groups smaller than this are sorted 
                                   // using insertion sort

// ----- "local" global variables
Int32 Shallow_limit;               // Max depth for shallow sorting
UChar *Shallow_text_limit;         // Text+Shallow_limit

#define UNROLL 1                   // if !=0 partially unroll shallow_mkq

// ----- some prototypes -------------
void helped_sort(Int32 *a, Int32 n, Int32 depth, Int32 limit);
static void shallow_inssort_lcp(Int32 *a, Int32 n, UChar *text_depth);

// --- static prototypes (gcc4 requires they are not inside functions) 
static void shallow_mkq(Int32 *a, int n, UChar *text_depth);
static void shallow_mkq16(Int32 *a, int n, UChar *text_depth);
static void shallow_mkq32(Int32 *a, int n, UChar *text_depth);


// ***** entry point for shallow sort routines *****
void shallow_sort(Int32 *a, int n, int shallow_limit) 
{ 
  // init global variables
  Shallow_limit = shallow_limit;        
  Shallow_text_limit = Text + shallow_limit;
  // call multikey quicksort
  // skip 2 chars since suffixes come from the same bucket 
  switch(_ds_Word_size) {
  case(1): shallow_mkq(a, n, Text+2); break;
  case(2): shallow_mkq16(a, n, Text+2); break;
  case(4): shallow_mkq32(a, n, Text+2); break;
  default:
    fprintf(stderr,
	    "Invalid word size for mkqs (%d) (shallow_sort)\n",_ds_Word_size);
    exit(1);
  }     
}


/* =======================================================
   auxiliary procedures and macro for bentley-sedgewick's
   multikey quicksort
   ======================================================= */
__inline__ void vecswap2(Int32 *a, Int32 *b, int n)
{   while (n-- > 0) {
        Int32 t = *a;
        *a++ = *b;
        *b++ = t;
    }
}

#define swap2(a, b) { t = *(a); *(a) = *(b); *(b) = t; }
#define ptr2char(i) (*(*(i) + text_depth))

__inline__ Int32 *med3func(Int32 *a, Int32 *b, Int32 *c, UChar *text_depth)
{   int va, vb, vc;
    if ((va=ptr2char(a)) == (vb=ptr2char(b)))
        return a;
    if ((vc=ptr2char(c)) == va || vc == vb)
        return c;       
    return va < vb ?
          (vb < vc ? b : (va < vc ? c : a ) )
        : (vb > vc ? b : (va < vc ? a : c ) );
}
#define med3(a, b, c) med3func(a, b, c, text_depth)


/* ********************************************************
   recursive multikey quicksort from Bentley-Sedgewick
   stops when text_depth reaches Shallow_depth_limit 
   that is when we have found that the current set of strings
   have Shallow_limit chars in common
   ******************************************************** */
static void shallow_mkq(Int32 *a, int n, UChar *text_depth)
{
  __inline__ void vecswap2(Int32 *a, Int32 *b, int n);
  int d, r, partval;
  Int32 *pa, *pb, *pc, *pd, *pl, *pm, *pn, t;
  UChar *next_depth;

  // ---- On small arrays use insertions sort
  if (n < Mk_qs_thresh) {
    shallow_inssort_lcp(a, n, text_depth);
    return;
  }

  // ----------- choose pivot --------------
 repeat:
  pl = a;
  pm = a + (n/2);
  pn = a + (n-1);
  if (n > 30) { // On big arrays, pseudomedian of 9
    d = (n/8);
    pl = med3(pl, pl+d, pl+2*d);
    pm = med3(pm-d, pm, pm+d);
    pn = med3(pn-2*d, pn-d, pn);
  }
  pm = med3(pl, pm, pn);
  swap2(a, pm);
  partval = ptr2char(a);
  pa = pb = a + 1;
  pc = pd = a + n-1;
  // -------- partition -----------------
  for (;;) {
    while (pb <= pc && (r = ptr2char(pb)-partval) <= 0) {
      if (r == 0) { swap2(pa, pb); pa++; }
      pb++;
    }
    while (pb <= pc && (r = ptr2char(pc)-partval) >= 0) {
      if (r == 0) { swap2(pc, pd); pd--; }
      pc--;
    }
    if (pb > pc) break;
    swap2(pb, pc);
    pb++;
    pc--;
  }

#if UNROLL
  if(pa>pd) {
    // all values were equal to partval: make it simpler
    if( (next_depth = text_depth+1) >= Shallow_text_limit) {
      helped_sort(a, n, next_depth-Text, Shallow_limit);
      return;
    }
    else {
      text_depth = next_depth;
      goto repeat;
    }
  }
#endif
  // partition a[] into the values smaller, equal, and larger that partval
  pn = a + n;
  r = min(pa-a, pb-pa);    vecswap2(a,  pb-r, r);
  r = min(pd-pc, pn-pd-1); vecswap2(pb, pn-r, r);
  // --- sort smaller strings -------
  if ((r = pb-pa) > 1)
    shallow_mkq(a, r, text_depth);
  // --- sort strings starting with partval -----
  if( (next_depth = text_depth+1) < Shallow_text_limit)
    shallow_mkq(a + r, pa-pd+n-1, next_depth);
  else 
    helped_sort(a + r, pa-pd+n-1, next_depth-Text, Shallow_limit);
  if ((r = pd-pc) > 1)
    shallow_mkq(a + n-r, r, text_depth);
}



/* ************** 16 *************** */
#define ptr2char16(i) (getword16(*(i) + text_depth))
#define getword16(s) ((unsigned)((*(s) << 8) | *((s)+1)))

#if 0
__inline__ Int32 *med3func16(Int32 *a, Int32 *b, Int32 *c, UChar *text_depth)
{   int va, vb, vc;
    if ((va=ptr2char16(a)) == (vb=ptr2char16(b)))
        return a;
    if ((vc=ptr2char16(c)) == va || vc == vb)
        return c;       
    return va < vb ?
          (vb < vc ? b : (va < vc ? c : a ) )
        : (vb > vc ? b : (va < vc ? a : c ) );
}
#define med3_16(a, b, c) med3func16(a, b, c, text_depth)
#endif

static void shallow_mkq16(Int32 *a, int n, UChar *text_depth)
{
  __inline__ void vecswap2(Int32 *a, Int32 *b, int n);
  int d, r, partval;
  Int32 *pa, *pb, *pc, *pd, *pl, *pm, *pn, t;
  UChar *next_depth;

  // ---- On small arrays use insertions sort
  if (n < Mk_qs_thresh) {
    shallow_inssort_lcp(a, n, text_depth);
    return;
  }

  // ----------- choose pivot --------------
 repeat:
  pl = a;
  pm = a + (n/2);
  pn = a + (n-1);
  if (n > 30) { // On big arrays, pseudomedian of 9
    d = (n/8);
    pl = med3(pl, pl+d, pl+2*d);
    pm = med3(pm-d, pm, pm+d);
    pn = med3(pn-2*d, pn-d, pn);
  }
  pm = med3(pl, pm, pn);
  swap2(a, pm);
  partval = ptr2char16(a);
  pa = pb = a + 1;
  pc = pd = a + n-1;
  // -------- partition -----------------
  for (;;) {
    while (pb <= pc && (r = ptr2char16(pb)-partval) <= 0) {
      if (r == 0) { swap2(pa, pb); pa++; }
      pb++;
    }
    while (pb <= pc && (r = ptr2char16(pc)-partval) >= 0) {
      if (r == 0) { swap2(pc, pd); pd--; }
      pc--;
    }
    if (pb > pc) break;
    swap2(pb, pc);
    pb++;
    pc--;
  }
#if UNROLL
  if(pa>pd) {
    // all values were equal to partval: make it simpler
    if( (next_depth = text_depth+2) >= Shallow_text_limit) {
      helped_sort(a, n, next_depth-Text, Shallow_limit);
      return;
    }
    else {
      text_depth = next_depth;
      goto repeat;
    }
  }
#endif
  // partition a[] into the values smaller, equal, and larger that partval
  pn = a + n;
  r = min(pa-a, pb-pa);    vecswap2(a,  pb-r, r);
  r = min(pd-pc, pn-pd-1); vecswap2(pb, pn-r, r);
  // --- sort smaller strings -------
  if ((r = pb-pa) > 1)
    shallow_mkq16(a, r, text_depth);
  // --- sort strings starting with partval -----
  if( (next_depth = text_depth+2) < Shallow_text_limit)
    shallow_mkq16(a + r, pa-pd+n-1, next_depth);
  else 
    helped_sort(a + r, pa-pd+n-1, next_depth-Text, Shallow_limit);
  if ((r = pd-pc) > 1)
    shallow_mkq16(a + n-r, r, text_depth);
}


/* *************** 32 **************** */
#define ptr2char32(i) (getword32(*(i) + text_depth))
#define getword32(s) ((unsigned)( (*(s) << 24) | ((*((s)+1)) << 16) \
                                  | ((*((s)+2)) << 8) | (*((s)+3)) ))
static void shallow_mkq32(Int32 *a, int n, UChar *text_depth)
{
  __inline__ void vecswap2(Int32 *a, Int32 *b, int n);
  UInt32 partval, val;
  Int32 *pa, *pb, *pc, *pd, *pl, *pm, *pn, t, d, r;
  UChar *next_depth;

  // ---- On small arrays use insertions sort
  if (n < Mk_qs_thresh) {
    shallow_inssort_lcp(a, n, text_depth);
    return;
  }

  // ----------- choose pivot --------------
 repeat:
  pl = a;
  pm = a + (n/2);
  pn = a + (n-1);
  if (n > 30) { // On big arrays, pseudomedian of 9
    d = (n/8);
    pl = med3(pl, pl+d, pl+2*d);
    pm = med3(pm-d, pm, pm+d);
    pn = med3(pn-2*d, pn-d, pn);
  }
  pm = med3(pl, pm, pn);
  swap2(a, pm);
  partval = ptr2char32(a);
  pa = pb = a + 1;
  pc = pd = a + n-1;
  // -------- partition -----------------
  for (;;) {
    while (pb <= pc &&  (val=ptr2char32(pb)) <=  partval) {
      if (val == partval) { swap2(pa, pb); pa++; }
      pb++;
    }
    while (pb <= pc && (val=ptr2char32(pc)) >= partval) {
      if (val == partval) { swap2(pc, pd); pd--; }
      pc--;
    }
    if (pb > pc) break;
    swap2(pb, pc);
    pb++;
    pc--;
  }
#if UNROLL
  if(pa>pd) {
    // all values were equal to partval: make it simpler
    if( (next_depth = text_depth+4) >= Shallow_text_limit) {
      helped_sort(a, n, next_depth-Text, Shallow_limit);
      return;
    }
    else {
      text_depth = next_depth;
      goto repeat;
    }
  }
#endif
  // partition a[] into the values smaller, equal, and larger that partval
  pn = a + n;
  r = min(pa-a, pb-pa);    vecswap2(a,  pb-r, r);
  r = min(pd-pc, pn-pd-1); vecswap2(pb, pn-r, r);
  // --- sort smaller strings -------
  if ((r = pb-pa) > 1)
    shallow_mkq32(a, r, text_depth);
  // --- sort strings starting with partval -----
  if( (next_depth = text_depth+4) < Shallow_text_limit)
    shallow_mkq32(a + r, pa-pd+n-1, next_depth);
  else 
    helped_sort(a + r, pa-pd+n-1, next_depth-Text, Shallow_limit);
  if ((r = pd-pc) > 1)
    shallow_mkq32(a + n-r, r, text_depth);
}



/* >>>>>>>>>>>>>>>>>>>>>> insertion sort routines >>>>>>>>>>>>>>>>>>>
   This insertion sort routines sorts the suffixes a[0] .. a[n-1]
   which have a common prexif of length text_depth-Text.
   The comparisons are done going at most at depth Shallow_limit;
   suffixes which have Shallow_limit chars in common are sorted using 
   helped_sort().
   This inserion_sort keeps trak of the lcp in order to speed up 
   the sorting.
  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */

/* ***********************************************************************
   Function to compare two strings originating from the *b1 and *b2
   The size of the unrolled loop must be at most equal to the costant 
   Cmp_overshoot defined in common.h
   When the function is called Cmp_left must contain the maximum number of 
   comparisons the algorithm can do before returning 0 (equal strings)
   At exit Cmp_left has been decreased by the # of comparisons done   
   *********************************************************************** */ 
static Int32 Cmp_left;
__inline__ 
Int32 cmp_unrolled_shallow_lcp(UChar *b1, UChar *b2)
{

  UChar c1, c2;
  assert(b1 != b2);

  // execute blocks of 16 comparisons until a difference
  // is found or we run out of the string 
  do {
    // 1
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 2
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_left -=  1; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 3
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_left -=  2; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 4
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_left -=  3; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 5
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_left -=  4; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 6
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_left -=  5; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 7
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_left -=  6; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 8
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_left -=  7; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 9
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_left -=  8; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 10
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_left -=  9; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 11
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_left -= 10; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 12
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_left -= 11; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 13
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_left -= 12; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 14
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_left -= 13; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 15
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_left -= 14; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 16
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_left -= 15; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // if we have done enough comparisons the strings are considered equal
    Cmp_left -= 16;
    if(Cmp_left<=0) return 0;
    // assert( b1<Upper_text_limit && b2<Upper_text_limit);
  } while(1);
  //return (b2-Text) - (b1-Text);   // we have  b2>b1 <=> *b2<*b1
  return b2 - b1;
} 

/* *****************************************************************
   this is the insertion sort routine called by multikey-quicksort
   for sorting small groups.
   During insertion sort the comparisons are done calling 
   cmp_unrolled_shallow_lcp() and two strings are equal if the coincides 
   for Shallow_limit characters.
   After this first phase we sort groups of "equal_string" using 
   helped_sort(). 
   Usage of lcp. 
   For i=1,...n-1 let lcp[i] denote the lcp between a[i] and a[i+1].
   assume a[0] ... a[j-1] are already ordered and that we want to 
   insert a new element ai. If suf(ai) >= suf(a[j-1]) we are done.
   If suf(ai)<suf(a[j-1]) we notice that: if lcpi>lcp[j-2] then
   suf(ai)>suf(a[j-2]) and we can stop since
   j-2 mmmmmmg
   j-1 mmmmmmmmmmmm
   ai  mmmmmmmmmmmf ] lcpi
   so we write a[j-1] in position j and ai in position j-1.
   if lcpi==lcp[j-2] then we need to compare suf(ai) with suf(a[j-2])
   j-2 mmmmmmmmmmm?     we can have either ?<f or ?>f or ?==f  
   j-1 mmmmmmmmmmmm
   j   mmmmmmmmmmmf
   so we move a[j-1] to position j and compare suf(ai) with suf(a[j-2])
   starting from lcpi.
   Finally, if lcpi<lcp[j-2] then
   j-2 mmmmmmmmmmmmmmmmmmg
   j-1 mmmmmmmmmmmmmmmmmmm
   j   mmmmmmmmmmmmmf
   hence we have suf(ai)<suf(a[j-2]) and we consider a[j-3];
   if lcpi<lcp[j-3] we go on look at a[j-4] and go on.
   if lcp[j]>lcp[j-3] we are in the following position:
   j-3 mmmmmmc                   
   j-2 mmmmmmmmmmmmmmmg        
   j-1 mmmmmmmmmmmmmmmm          
   j   mmmmmmmmmmf
   and we know that suf(ai) is larger than suf(a[j-3]). If we find that 
   lcpi==lcp[j-3] then we must compare suf(ai) with suf(a[j-3])
   but starting with position lcpi
   ***************************************************************** */
static int lcp_aux[1+Max_thresh];
static int *lcp=lcp_aux+1; 
static void shallow_inssort_lcp(Int32 *a, Int32 n, UChar *text_depth)
{   
  __inline__ Int32 cmp_unrolled_shallow_lcp(UChar *, UChar *);
  Int32 i, j, j1, lcp_new, r, ai,lcpi;
  Int32 cmp_from_limit;
  UChar *text_depth_ai;

  // --------- initialize ----------------
  lcp_aux[0] = -1;               // set lcp[-1] = -1
  for(i=0;i<n;i++) lcp[i]=0;     // I think this loop is not necessary
  // cmp_from_limit is # of cmp's to be done to reach Shallow_limit cmp's
  cmp_from_limit = Shallow_text_limit-text_depth;

  // ----- start insertion sort -----------
  for (i = 1; i< n ; i++) {
    ai = a[i]; lcpi = 0;
    text_depth_ai = ai + text_depth;
    j=i; j1=j-1;                  // j1 is a shorhand for j-1
    while(1) {           

      // ------ compare ai with a[j-1] --------
      Cmp_left = cmp_from_limit-lcpi;  
      r = cmp_unrolled_shallow_lcp(lcpi+a[j1]+text_depth,lcpi+text_depth_ai);
      lcp_new = cmp_from_limit - Cmp_left;       // lcp between ai and a[j1] 
      assert(r!=0 || lcp_new>= cmp_from_limit);

      if(r<=0) {         // we have a[j-1] <= ai
	lcp[j1]=lcp_new; // ai will be written in a[j]; update lcp[j-1]
	break;
      }

      // --- we have a[j-1]>ai. a[j-1] and maybe other will be moved down 
      // --- use lcp to move down as many elements of a[] as possible
      lcpi = lcp_new;                
      do {
	a[j] = a[j1];               // move down a[j-1]
        lcp[j] = lcp[j1];           // move down lcp[j-1]
	j=j1; j1--;                 // update j and j1=j-1
      } while(lcpi<lcp[j1]);        // recall that lcp[-1]=-1

      if(lcpi>lcp[j1]) break;       // ai will be written in position j

      // if we get here lcpi==lcp[j1]: we will compare them at next iteration

    }     // end for(j=i ...
    a[j]=ai;
    lcp[j]=lcpi;
  }       // end for(i=1 ... 

  // ----- done with insertion sort. now sort groups of equal strings
  for(i=0;i<n-1;i=j+1) {
    for(j=i; j<n ;j++)
      if(lcp[j]<cmp_from_limit) break;
    if(j-i>0) {
      helped_sort(a+i,j-i+1,Shallow_limit, Shallow_limit); 
    }
  }
}











