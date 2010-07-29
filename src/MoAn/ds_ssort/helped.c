/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   helped.c  
   these routines sort a group of strings which have a common prefix
   using induced sorting (if possible) or a deep-sort routine
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include "common.h"

// ----------------- external variables ---------------------
extern Int32 Anchor_dist;       // Distance between anchors
extern Int32 Anchor_num;        // Number of anchor points
extern Int32 *Anchor_rank;      // rank (in the sorted suffixes of the  
                                // anchor points (-1 if rank is unknown)
extern UInt16 *Anchor_offset;   // offset (wrt to the anchor) of the suffix
                                // whose rank is in Anchor_rank. 

extern UChar *Text;               // start of input string
extern Int32 ftab[];              // table of buckets endpoints
extern Int32 *Sa;                 // suffix array
extern Int32 Max_pseudo_anchor_offset; // maximum offset considered when 
                                       // searching a pseudo anchor
extern Int32 B2g_ratio;           // maximum ratio bucket_size/group_size
                                  // accepted for pseudo anchor_sorting
extern Int32 Update_anchor_ranks; // if!=0 update anchor ranks when determining
                                  // rank for pseudo-sorting

extern Int32 Calls_helped_sort;     
extern Int32 Calls_anchor_sort_forw;     
extern Int32 Calls_anchor_sort_backw;     
extern Int32 Calls_pseudo_anchor_sort_forw;     

// function for deep sorting the suffixes in a[0] .. a[n-1]
// deep is the length of the prefixes known to be equal in all suffixes
void deep_sort(Int32 *a, Int32 n, Int32 depth);

// macro to compute the bucket for the suffix
// starting at pos. Note that since pos is evaluated twice
// it should be an expression without side-effects
#define Get_small_bucket(pos) ((Text[pos]<<8) + Text[pos+1])

// --- static prototypes (gcc4 requires they are not inside functions) 
static 
void general_anchor_sort(Int32 *a,Int32 n,Int32 pos,Int32 rank,Int32 off);
static 
void pseudo_anchor_sort(Int32 *a, Int32 n,Int32 pseudo_an,Int32 offset);
static Int32 split_group(Int32 *a, int n, int,int,Int32,int *);
static void update_anchors(Int32 *a, Int32 n);
static void pseudo_or_deep_sort(Int32 *a, Int32 n, Int32 depth);



/* *****************************************************************
   This procedure sort the strings a[0] ... a[n-1] with the help of an
   anchor. The real sorting is done by the procedure
   anchor_sort(). Here we choose the anchor.  The parameter depth is
   the number of chars that a[0] ... a[n-1] are known to have in
   common (thus a direct comparison among a[i] and a[j] should start
   from position depth) Note that a[] is a subsection of the sa therefore
   a[0] ... a[n-1] are starting position of suffixes
   For every a[i] we look at the anchor a[i]/Anchor_dist and the one 
   after that. This justifies the definition of Anchor_num (the size of
   Anchor_ofset[] and Anchor_rank[] defined in ds_sort()) as
     Anchor_num = 2 + (n-1)/Anchor_dist    
   ***************************************************************** */
void helped_sort(Int32 *a, int n, int depth, int limit)
{ 
  Int32 i, curr_sb, diff, toffset, aoffset;
  Int32 text_pos, anchor_pos, anchor, anchor_rank;
  Int32 min_forw_offset, min_forw_offset_buc, max_back_offset;
  Int32 best_forw_anchor, best_forw_anchor_buc, best_back_anchor; 
  Int32 forw_anchor_index, forw_anchor_index_buc, back_anchor_index;


  if(depth >= limit)
  {
  	 
     return;
  }
 printf("Helped invoked with depth < limit: %d", depth);
 Calls_helped_sort++;          // update count

  
  if(n==1) goto done_sorting;    // simplest case: only one string

  // if there are no anchors use pseudo-anchors or deep_sort
  if(Anchor_dist==0) {
    pseudo_or_deep_sort(a, n, depth);
    return;
  }

  // compute the current bucket 
  curr_sb = Get_small_bucket(a[0]);

  // init best anchor variables with illegal values
  min_forw_offset = min_forw_offset_buc = INT_MAX;
  max_back_offset = INT_MIN;
  best_forw_anchor = best_forw_anchor_buc = best_back_anchor = -1; 
  forw_anchor_index = forw_anchor_index_buc = back_anchor_index = -1;
  // look at the anchor preceeding each a[i]
  for(i=0;i<n;i++) {
    text_pos = a[i];
    // get anchor preceeding text_pos=a[i]
    anchor = text_pos/Anchor_dist;
    toffset = text_pos % Anchor_dist;  // distance of a[i] from anchor
    aoffset = Anchor_offset[anchor];   // distance of sorted suf from anchor 
    if(aoffset<Anchor_dist) {          // check if it is a "sorted" anchor
      diff = aoffset - toffset;
      assert(diff!=0);
      if(diff>0) {     // anchor <=  a[i] < (sorted suffix)
	if(curr_sb!=Get_small_bucket(text_pos+diff)) {
	  if(diff<min_forw_offset) {
	    min_forw_offset = diff;
	    best_forw_anchor = anchor;
	    forw_anchor_index = i;
	  }
	}
	else {  // the sorted suffix belongs to the same bucket of a[0]..a[n-1]
	  if(diff<min_forw_offset_buc) {
	    min_forw_offset_buc = diff;
	    best_forw_anchor_buc = anchor;
	    forw_anchor_index_buc = i;
	  }
	}
      }
      else {          // diff<0 =>  anchor <= (sorted suffix) < a[i]
	if(diff>max_back_offset) {
	  max_back_offset = diff;
	  best_back_anchor = anchor;
	  back_anchor_index = i;
	}
	// try to find a sorted suffix > a[i] by looking at next anchor
	aoffset = Anchor_offset[++anchor];
	if(aoffset<Anchor_dist) {
	  diff = Anchor_dist + aoffset - toffset;
	  assert(diff>0);
	  if(curr_sb!=Get_small_bucket(text_pos+diff)) {
	    if(diff<min_forw_offset) {
	      min_forw_offset = diff;
	      best_forw_anchor = anchor;
	      forw_anchor_index = i;
	    }
	  } else {
	    if(diff<min_forw_offset_buc) {
	      min_forw_offset_buc = diff;
	      best_forw_anchor_buc = anchor;
	      forw_anchor_index_buc = i;
	    }
	  }
	}
      }
    }
  }
  // ------ if forward anchor_sort is possible, do it! --------	    
  if(best_forw_anchor>=0 && min_forw_offset<depth-1) {
    Calls_anchor_sort_forw++;
    assert(min_forw_offset<2*Anchor_dist);
    anchor_pos = a[forw_anchor_index] + min_forw_offset;
    anchor_rank = Anchor_rank[best_forw_anchor];
    assert(Sa[anchor_rank]==anchor_pos);
    general_anchor_sort(a,n,anchor_pos,anchor_rank,min_forw_offset);
    goto done_sorting;
  }
  // ------ if backward anchor_sort is possible do it! ---------
  if(best_back_anchor>=0) {
    UChar *T0, *Ti; int j;

    assert(max_back_offset>-Anchor_dist && max_back_offset<0);
    // make sure that the offset is legal for all a[i]
    for(i=0;i<n;i++) {
      if(a[i]+max_back_offset<0) 
	goto fail;                    // illegal offset, give up
    }
    // make sure that a[0] .. a[n-1] are preceded by the same substring
    T0 = Text + a[0];
    for(i=1;i<n;i++) {
      Ti = Text + a[i];
      for(j=max_back_offset; j<= -1; j++)
	if(T0[j]!=Ti[j]) goto fail;   // mismatch, give up
    }
    // backward anchor sorting is possible
    Calls_anchor_sort_backw++;
    anchor_pos = a[back_anchor_index] + max_back_offset;
    anchor_rank = Anchor_rank[best_back_anchor];
    assert(Sa[anchor_rank]==anchor_pos);
    general_anchor_sort(a,n,anchor_pos,anchor_rank,max_back_offset);
    goto done_sorting;
  }
 fail:
  // ----- try forward anchor_sort with anchor in the same bucket
  if(best_forw_anchor_buc>=0 && min_forw_offset_buc<depth-1) {
    int equal,lower,upper;

    assert(min_forw_offset_buc<2*Anchor_dist);
    anchor_pos = a[forw_anchor_index_buc] + min_forw_offset_buc;
    anchor_rank = Anchor_rank[best_forw_anchor_buc];
    assert(Sa[anchor_rank]==anchor_pos);

    // establish how many suffixes can be sorted using anchor_sort()
    equal=split_group(a,n,depth,min_forw_offset_buc,
                                forw_anchor_index_buc,&lower);
    if(equal==n) {
      Calls_anchor_sort_forw++;
      general_anchor_sort(a,n,anchor_pos,anchor_rank,min_forw_offset_buc);
    }
    else {
      //  -- a[0] ... a[n-1] are split into 3 groups: lower, equal, upper
      upper = n-equal-lower;
      assert(upper>=0);
      // printf("Warning! lo=%d eq=%d up=%d a=%x\n",lower,equal,upper,(int)a);
      // sort the equal group 
      Calls_anchor_sort_forw++;
      if(equal>1)
	general_anchor_sort(a+lower,equal,anchor_pos,anchor_rank,
			    min_forw_offset_buc);

      // sort upper and lower groups using deep_sort
      if(lower>1) pseudo_or_deep_sort(a,lower,depth);
      if(upper>1) pseudo_or_deep_sort(a+lower+equal,upper,depth);
    }       // end if(equal==n) ... else
    goto done_sorting;
  }         // end hard case

  // ---------------------------------------------------------------
  // If we get here it means that everything failed
  // In this case we simply deep_sort a[0] ... a[n-1]
  // ---------------------------------------------------------------
  pseudo_or_deep_sort(a, n, depth);
 done_sorting:
  // -------- update Anchor_rank[], Anchor_offset[] ------- 
  if(Anchor_dist>0) update_anchors(a, n);
}
  


/* *******************************************************************
   try pseudo_anchor sort or deep_sort
   ******************************************************************** */
static void pseudo_or_deep_sort(Int32 *a, Int32 n, Int32 depth)
{

  Int32 offset, text_pos, sb, pseudo_anchor_pos, max_offset, size;
 
  // ------- search for a useful pseudo-anchor -------------
  if(Max_pseudo_anchor_offset>0) {

    max_offset = min(depth-1,Max_pseudo_anchor_offset);
    text_pos = a[0];
    for(offset=1;offset<max_offset;offset++) {
      pseudo_anchor_pos = text_pos+offset;
      sb = Get_small_bucket(pseudo_anchor_pos);
      // check if pseudo_anchor is in a sorted bucket
      if(IS_SORTED_BUCKET(sb)) {
	size=BUCKET_SIZE(sb);                     // size of group
	if(size>B2g_ratio*n) continue;            // discard large groups 
	// sort a[0] ... a[n-1] using pseudo_anchor
	pseudo_anchor_sort(a,n,pseudo_anchor_pos,offset);
	Calls_pseudo_anchor_sort_forw++;        // update count
	return;
      }
    }
  }
  deep_sort(a,n,depth);
}

/* ********************************************************************
   this routine sorts the suffixes a[0] ... a[n-1] using the fact that
   in their common prefix, after offset characters, there is a 
   suffix which is in an already sorted bucket. This suffix is called
   a pseudo anchor since it is used essentially as an anchor, but
   it is not in an anchor position (=position multiple of Anchor_dist)
   ******************************************************************** */
static 
void pseudo_anchor_sort(Int32 *a,Int32 n,Int32 pseudo_anchor_pos, Int32 offset)
{
  Int32 get_rank(Int32);
  Int32 get_rank_update_anchors(Int32);
  Int32 pseudo_anchor_rank;

  // ---------- compute rank ------------
  if(Update_anchor_ranks!=0 && Anchor_dist>0)
    pseudo_anchor_rank = get_rank_update_anchors(pseudo_anchor_pos);
  else
    pseudo_anchor_rank = get_rank(pseudo_anchor_pos);
  // ---------- check rank --------------
  assert(Sa[pseudo_anchor_rank]==pseudo_anchor_pos);
  // ---------- do the sorting ----------
  general_anchor_sort(a,n,pseudo_anchor_pos,pseudo_anchor_rank,offset);
}


/* ********************************************************
   macros for marking integers: works assuming integers have 
   at least 32 bit and that the 32nd bit is not used
   This simply means that the text size can be at most 2GB
   ********************************************************* */
#define MARKER (1<<31)
#define MARK(i) {                \
  assert(( Sa[i]&MARKER) == 0);  \
  (Sa[i] |= MARKER);             \
}
#define ISMARKED(i) (Sa[i] & MARKER)
#define UNMARK(i) (Sa[i] &= ~MARKER)

/* ********************************************************************
   This routines sorts a[0] ... a[n-1] using the fact that
   in their common prefix, after offset characters, there is a 
   suffix whose rank is known. In this routine we call this suffix anchor
   (and we denote its position and rank with anchor_pos and anchor_rank 
   respectively) but it is not necessarily an anchor (=does not necessarily 
   starts at position multiple of Anchor_dist) since this function is
   called by pseudo_anchor_sort().
   The routine works by scanning the suffixes before and after the anchor
   in order to find (and mark) those which are suffixes of a[0] ... a[n-1].
   After that, the ordering of a[0] ... a[n-1] is derived with a sigle
   scan of the marked suffixes.
   ******************************************************************** */
static void general_anchor_sort(Int32 *a, Int32 n, 
                         Int32 anchor_pos, Int32 anchor_rank, Int32 offset)
{
  int integer_cmp(const void *, const void *);
  Int32 sb, lo, hi;
  Int32 curr_lo, curr_hi, to_be_found, i,j;
  Int32 item; 
  void *ris;

  assert(Sa[anchor_rank]==anchor_pos);
  /* ---------- get bucket of anchor ---------- */
  sb = Get_small_bucket(anchor_pos);
  lo = BUCKET_FIRST(sb);
  hi = BUCKET_LAST(sb);
  assert(sb==Get_small_bucket(a[0]+offset));
  // ------ sort pointers a[0] ... a[n-1] as plain integers
  qsort(a,n, sizeof(Int32), integer_cmp);

  // ------------------------------------------------------------------
  // now we scan the bucket containing the anchor in search of suffixes
  // corresponding to the ones we have to sort. When we find one of
  // such suffixes we mark it. We go on untill n sfx's have been marked 
  // ------------------------------------------------------------------
  curr_hi = curr_lo = anchor_rank;

  // the anchor must correspond to a suffix to be sorted
  #if DEBUG
  item = anchor_pos-offset;
  assert(bsearch(&item,a,n,sizeof(Int32), integer_cmp));
  #endif

  MARK(curr_lo);
  // scan suffixes preceeding and following the anchor
  for(to_be_found=n-1;to_be_found>0; ) {
    // invariant: the next positions to check are curr_lo-1 and curr_hi+1
    assert(curr_lo > lo || curr_hi < hi);
    while (curr_lo > lo) {
      item = Sa[--curr_lo]-offset;
      ris = bsearch(&item,a,n,sizeof(Int32), integer_cmp);
      if(ris)	{MARK(curr_lo); to_be_found--;}
      else	break;
    }
    while (curr_hi < hi) {
      item = Sa[++curr_hi]-offset;
      ris = bsearch(&item,a,n,sizeof(Int32), integer_cmp);
      if(ris)	{MARK(curr_hi); to_be_found--;}
      else      break;
    }
  }
  // sort a[] using the marked suffixes
  for(j=0, i=curr_lo;i<=curr_hi;i++) 
    if(ISMARKED(i)) {
      UNMARK(i);
      a[j++] = Sa[i] - offset;
    }
  assert(j==n);  // make sure n items have been sorted
}

/* ********************************************************************
   compute the rank of the suffix starting at pos.
   It is required that the suffix is in an already sorted bucket
   ******************************************************************** */
Int32 get_rank(Int32 pos)
{
  Int32 sb, lo, hi, j;

  sb = Get_small_bucket(pos);  
  if(!IS_SORTED_BUCKET(sb)) {
    fprintf(stderr,"Illegal call to get_rank! (get_rank1)\n");
    exit(1);
  }
  lo = BUCKET_FIRST(sb);
  hi = BUCKET_LAST(sb);
  for(j=lo;j<=hi;j++) 
    if(Sa[j]==pos) return j;
  fprintf(stderr,"Illegal call to get_rank! (get_rank2)\n");
  exit(1);
  return 1;   // so that the compiler does not complain
}

/* ********************************************************************
   compute the rank of the suffix starting at pos. At the same time
   check if the rank of the suffixes in the bucket containing pos
   can be used to update some entries in Anchor_offset[] and Anchor_rank[]
   It is required that the suffix is in an already sorted bucket   
   ******************************************************************** */
static UChar bucket_ranked[65536];
Int32 get_rank_update_anchors(Int32 pos)
{
  Int32 get_rank(Int32 pos);
  Int32 sb, lo, hi, j, toffset, aoffset, anchor, rank;

  assert(Anchor_dist>0);
  // --- get bucket and verify it is a sorted one
  sb = Get_small_bucket(pos);  
  if(!(IS_SORTED_BUCKET(sb))) {
    fprintf(stderr,"Illegal call to get_rank! (get_rank_update_anchors)\n");
    exit(1);
  }
  // --- if the bucket has been already ranked just compute rank; 
  if(bucket_ranked[sb]) return get_rank(pos);
  // --- rank all the bucket 
  bucket_ranked[sb]=1;
  rank = -1;
  lo = BUCKET_FIRST(sb);
  hi = BUCKET_LAST(sb);
  for(j=lo;j<=hi;j++) {  
    // see if we can update an anchor
    toffset = Sa[j]%Anchor_dist;
    anchor  = Sa[j]/Anchor_dist;
    aoffset = Anchor_offset[anchor];  // dist of sorted suf from anchor 
    if(toffset<aoffset) {
      Anchor_offset[anchor] = toffset;
      Anchor_rank[anchor] = j;
    }
    // see if we have found the rank of pos, if so store it in rank
    if(Sa[j]==pos) {
      assert(rank==-1); rank=j;
    }
  }
  assert(rank>=0);
  return rank;
}


/* ****************************************************************** 
   comparison function used to sort pointers as if they were integers
   if a pointer does not correspond to an Int32 we must change the
   function accordingly 
   ****************************************************************** */
int integer_cmp(const void *a, const void *b)
{
  return *((Int32 *) a) -  *((Int32 *) b); 
}

/* ****************************************************************
   given a SORTED array of suffixes a[0] .. a[n-1]
   updates Anchor_rank[] and Anchor_offset[]
   **************************************************************** */
static void update_anchors(Int32 *a, Int32 n)
{
  Int32 i,anchor,toffset,aoffset,text_pos;

  assert(Anchor_dist>0);
  for(i=0;i<n;i++) {
    text_pos = a[i];
    // get anchor preceeding text_pos=a[i]
    anchor = text_pos/Anchor_dist;
    toffset = text_pos % Anchor_dist;     // distance of a[i] from anchor
    aoffset = Anchor_offset[anchor];  // dist of sorted suf from anchor 
    if(toffset<aoffset) {
      Anchor_offset[anchor] = toffset;
      Anchor_rank[anchor] = (a - Sa) + i;
      assert(Sa[Anchor_rank[anchor]]==
	     anchor*Anchor_dist+Anchor_offset[anchor]);
    }
  }
}
 


/* *******************************************************************
   This function takes as input an array a[0] .. a[n-1] of suffixes
   which share the first "depth" chars. "pivot" in an index in 0..n-1
   and offset and integer>0. The function splits a[0] .. a[n-1]
   into 3 groups: first the suffixes which are smaller than a[pivot],
   then those which are equal to a[pivot] and finally those which are 
   greater than a[pivot]. Here, smaller, equal, larger refer to 
   a lexicographic ordering limited to the first depth+offest chars
   (since the first depth chars are equal we only look at the chars
   in position depth, depth+1, ... depth+offset-1).
   The function returns the number "num" of suffixes equal to a[pivot],
   and stores in *first the first of these suffixes. So at the end 
   the smaller suffixes are in a[0] ... a[first-1],
   the equal suffixes in a[first] ... a[first+num-1],
   the larger suffixes in a[first+num] ... a[n-1]
   The splitting is done using a modified mkq()
   ******************************************************************* */
#define swap2(a, b) { t = *(a); *(a) = *(b); *(b) = t; }
#define ptr2char(i) (*(*(i) + text_depth))
static 
Int32 split_group(Int32 *a, int n, int depth,int offset,Int32 pivot,int *first)
{
  void vecswap2(Int32 *a, Int32 *b, int n);
  int r, partval;
  Int32 *pa, *pb, *pc, *pd, *pa_old, *pd_old, pivot_pos, t;
  UChar *text_depth,*text_limit;

  // --------- initialization ------------------------------------
  pivot_pos = a[pivot];       // starting position in T[] of pivot
  text_depth = Text+depth;
  text_limit = text_depth+offset;

  // -------------------------------------------------------------
  // In the following for() loop:
  // [pa ... pd] is the current working region, 
  // pb moves from pa towards pd 
  // pc moves from pd towards pa
  // -------------------------------------------------------------
  pa = a; pd = a + n-1;

  for(  ; pa!=pd && (text_depth<text_limit); text_depth++) {
    assert(pa<pd);
    // ------ the pivot char is Text[pivot_pos+depth] where 
    // depth = text_depth-Text. This is text_depth[pivot_pos]
    partval = text_depth[pivot_pos];
    // ----- partition ------------ 
    pb = pa_old = pa; 
    pc = pd_old = pd; 
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
    r = min(pa-pa_old, pb-pa); vecswap2(pa_old,  pb-r, r);
    r = min(pd-pc, pd_old-pd); vecswap2(pb, pd_old+1-r, r);
    // ------ compute new boundaries ----- 
    pa = pa_old + (pb-pa);     // there are pb-pa chars < partval
    pd = pd_old - (pd-pc);     // there are pd-pc chars > partval

  }
  *first=pa-a;        // index in a[] of the first suf. equal to pivot
  assert(pd-pa>=0);
  return pd-pa+1;     // return number of suffixes equal to pivot
}







