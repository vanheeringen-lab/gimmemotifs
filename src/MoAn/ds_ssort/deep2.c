/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   deep.c

   "deep" sorting routines
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "common.h"


/* ------ external global variables ------- */
extern UChar  *Text;                   // input string+ overshoot
extern UChar  *Upper_text_limit;       // Text+Text_size
extern Int32  Text_size;               // size of input string
extern Int32  Blind_sort_ratio;        // ratio for using blind_sort
extern Int32 Calls_deep_sort;     

/* ***********************************************************************
   Function to compare two strings originating from the *b1 and *b2
   The size of the unrolled loop must be at most equal to the costant 
   Cmp_overshoot defined in common.h
   the function return the result of the comparison (+ or -) and writes 
   in Cmp_done the number of successfull comparisons done
   *********************************************************************** */ 
static Int32 Cmp_done;
__inline__
Int32 cmp_unrolled_lcp(UChar *b1, UChar *b2)
{

  UChar c1, c2;
  assert(b1 != b2);
  Cmp_done=0;

  // execute blocks of 16 comparisons untill a difference
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
      Cmp_done +=  1; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 3
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_done +=  2; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 4
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_done +=  3; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 5
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_done +=  4; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 6
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_done +=  5; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 7
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_done +=  6; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 8
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_done +=  7; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 9
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_done +=  8; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 10
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_done +=  9; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 11
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_done += 10; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 12
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_done += 11; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 13
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_done += 12; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 14
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_done += 13; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 15
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_done += 14; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 16
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_done += 15; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 

    Cmp_done += 16;

  } while( b1<Upper_text_limit && b2<Upper_text_limit);

  //return (b2-Text) - (b1-Text);   // we have  b2>b1 <=> *b2<*b1
  return b2 - b1;
} 

/* **************************************************************
   ternary quicksort (seward-like) with lcp information
   ************************************************************** */
#define STACK_SIZE 100
#define Swap(i,j) {tmp=a[i]; a[i]=a[j]; a[j]=tmp;}
#define Pushd(x,y,z) {stack_lo[sp]=x; stack_hi[sp]=y; stack_d[sp]=z; sp++;}
#define Popd(x,y,z)  {sp--; x=stack_lo[sp]; y=stack_hi[sp]; z=stack_d[sp];} 
void qs_unrolled_lcp(Int32 *a, int n, int depth, int blind_limit)
{ 
  void blind_ssort(Int32 *a, Int32 n, Int32 depth);
  __inline__ Int32 cmp_unrolled_lcp( UChar *b1, UChar *b2 );
  UChar *text_depth, *text_pos_pivot;
  Int32 stack_lo[STACK_SIZE];
  Int32 stack_hi[STACK_SIZE];
  Int32 stack_d[STACK_SIZE];
  Int32 sp,r,r3,med,tmp;
  Int32 i, j, lo, hi,ris,lcp_lo,lcp_hi;

  // ----- init quicksort --------------
  r=sp=0;
  Pushd(0,n-1,depth);

  // ----- repeat untill stack is empty ------
  while (sp > 0) {
    assert ( sp < STACK_SIZE );
    Popd(lo,hi,depth);
    text_depth = Text+depth;

    // --- use shellsort for small groups
    if(hi-lo<blind_limit) { 
       blind_ssort(a+lo,hi-lo+1,depth);
       continue;
    }

    /* Random partitioning. Guidance for the magic constants 
       7621 and 32768 is taken from Sedgewick's algorithms
       book, chapter 35.
    */
    r = ((r * 7621) + 1) % 32768;
    r3 = r % 3;
    if (r3 == 0) med = lo; else
    if (r3 == 1) med = (lo+hi)>>1; else
                 med = hi;

    // --- partition ----
    Swap(med,hi);  // put the pivot at the right-end
    text_pos_pivot=text_depth+a[hi];
    i=lo-1; j=hi; 
    lcp_lo=lcp_hi=INT_MAX;
    while(1) {
      while(++i<hi) {
	ris=cmp_unrolled_lcp(text_depth+a[i], text_pos_pivot);
        if(ris>0) {
	  if(Cmp_done < lcp_hi) lcp_hi=Cmp_done; break;
	} else if(Cmp_done < lcp_lo) lcp_lo=Cmp_done;
      }
      while(--j>lo) {
	ris=cmp_unrolled_lcp(text_depth+a[j], text_pos_pivot);
        if(ris<0) { if(Cmp_done < lcp_lo) lcp_lo=Cmp_done; break; }
	else if(Cmp_done < lcp_hi) lcp_hi=Cmp_done;
      }
      if (i >= j) break; 
      Swap(i,j);
    }
    Swap(i,hi);  // put pivot at the middle

    // ---- testing ---------
    assert(lcp_lo<INT_MAX || i==lo);
    assert(lcp_hi<INT_MAX || i==hi);

    // --------- insert subproblems in stack; smallest last
    if(i-lo < hi-i) {
      Pushd(i+1,hi,depth+lcp_hi);
      if(i-lo>1) Pushd(lo,i-1,depth+lcp_lo);
    }
    else {
      Pushd(lo,i-1,depth+lcp_lo);
      if(hi-i>1) Pushd(i+1,hi,depth+lcp_hi);
    }
  }
}



/* ****************************************************************
   routine for deep-sorting the suffixes a[0] ... a[n-1]
   knowing that they have a common prefix of length "depth"
  **************************************************************** */   
void deep_sort(Int32 *a, Int32 n, Int32 depth)
{
  void blind_ssort(Int32 *a, Int32 n, Int32 depth);
  int blind_limit;

  Calls_deep_sort++;    
  assert(n>1);    // test to discover useless calls

  blind_limit=Text_size/Blind_sort_ratio;
  if(n<=blind_limit)
    blind_ssort(a,n,depth);  // small_group
  else 
    qs_unrolled_lcp(a,n,depth,blind_limit);
}














