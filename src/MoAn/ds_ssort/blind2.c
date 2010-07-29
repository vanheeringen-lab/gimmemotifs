/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   string sorting routine based on a blind tree
   26-jun-01 ver 1.0
   03-jul-01 ver 1.1 (get rid of node header)
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */ 
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "common.h"

/* ==================================================================
   comment: it is a little triky how we handle the case in which 
   there are two strings s1 s2 such that s1 is a prefix of s2.
   (the correct ordering is that s1 preceeds lexicographically s2)
   We proceed as follows. We insert the strings in order of increasing 
   length so that s1 is inserted before s2. When s2 is inserted we put
   it in a leaf which is to the left of s1's leaf (obviously they have
   the same parent node). This is wrong acording to the alphabetic
   ordering but is done so that is there is a third string s3 which
   has s1 as a prefix we are certain that s3 meets s2's leaf and not s1's.
   When we traverse the trie to get the sorted string we check if there 
   are two sibling with the same key and if so we invert them 
   to get the correct ordering
   =================================================================== */

/* ------ external global variables ------- */
extern UChar  *Text;                   // input string+ overshoot
extern Int32  Text_size;               // size of input string
extern UChar  *Upper_text_limit;       // Text+Text_size

/* ------- node of blind trie -------- */ 
typedef struct nodex {
  Int32 skip;
  UChar key;  
  struct nodex  *down;      // first child
  struct nodex *right;      // next brother
} node;


// -------- local global variables -------------------
#define BUFSIZE 1000
#define FREESIZE 5000
void *freearr[FREESIZE];
static node *bufn;
static int bufn_num=0, free_num=0;
static Int32 *Aux, Aux_written;
node **Stack;
int Stack_size;

// --- static prototypes (gcc4 requires they are not inside functions) 
static int neg_integer_cmp(const void *, const void *);
static node *find_companion(node *head, UChar *s);
static void insert_suffix(node *h, Int32 suf, int n, UChar mmchar);
static void traverse_trie(node *h);
static Int32 compare_suffixes(Int32 suf1, Int32 suf2, Int32 depth);
static void free_node_mem();
static node *get_leaf(node *head);



/* ****************************************************************
   routine for deep-sorting the suffixes a[0] ... a[n-1]
   knowing that they have a common prefix of length "depth"
  **************************************************************** */   
void blind_ssort(Int32 *a, Int32 n, Int32 depth)
{
  Int32 i,j,aj,lcp;
  node nh, *root, *h;

  // ---- sort suffixes in order of increasing length
  qsort(a,n, sizeof(Int32), neg_integer_cmp);

  // --- skip suffixes which have already reached the end-of-text
  for(j=0;j<n;j++)
    if(a[j]+depth < Text_size)
      break;
  if(j>=n-1) return;  // everything is already sorted!

  // ------ init stack -------
  Stack = (node **) malloc(n*sizeof(node *));
  if(Stack==NULL) {
    fprintf(stderr,"Out of memory! (blind_ssort)\n");
    exit(1);
  }

  // ------- init root with the first unsorted suffix
  nh.skip = -1;   nh.right = NULL; nh.down = (void *) a[j]; 
  root = &nh;

  // ------- insert suffixes a[j+1] ... a[n-1]
  for(i=j+1;i<n;i++) {
    h=find_companion(root, Text+a[i]);
    assert(h->skip==-1);
    assert(Stack_size<=i-j);
    aj=(Int32) h->down;
    assert(aj>a[i]);
    lcp = compare_suffixes(aj,a[i],depth);
    insert_suffix(root, a[i], lcp, Text[aj+lcp]);
  }

  // ---- traverse the trie and get suffixes in lexicographic order  
  Aux=a;  Aux_written = j;
  traverse_trie(root);
  assert(Aux_written==n);
 
  free_node_mem();
  free(Stack);
}

/* ***********************************************************************
   this function traverses the trie rooted at head following the string s. 
   Returns the leaf "corresponding" to the string s
   *********************************************************************** */
static node *find_companion(node *head, UChar *s)
{
  UChar c;
  node *p;
  int t;

  Stack_size = 0;                // init stack
  while(head->skip >= 0) {
    Stack[Stack_size++] = head;
    t = head->skip;
    if(s+t>=Upper_text_limit)    // s[t] does not exist: mismatch 
      return get_leaf(head);
    c = s[t]; p = head->down;
  repeat:
    if(c==p->key) {              // found branch corresponding to c
      head = p;
      continue;
    }
    else if(c<p->key)            // no branch corresponding to c: mismatch
      return get_leaf(head);
    if((p=(p->right))==NULL)     // no other branches: mismatch
      return get_leaf(head);
    goto repeat;                 // look at next branch
  }
  Stack[Stack_size++] = head;
  return head;
}


// this function returns a leaf below "head". 
// any leaf will do for the algorithm: we take the easiest to reach
static node *get_leaf(node *head)
{
  assert(head->skip>=0);

  do {
    head = head->down;
  } while(head->skip>=0);
  return head;
}



__inline__ node *new_node__blind_ssort(void)
{
  if(bufn_num-- == 0) {
    bufn = (node *) malloc(BUFSIZE * sizeof(node));
    if(bufn==NULL) {
      fprintf(stderr,"Out of mem (new_node1)\n"); exit(1);}
    freearr[free_num++] = (void *) bufn; 
    if(free_num>=FREESIZE) {
      fprintf(stderr,"Out of mem (new_node2)\n"); exit(1);}
   bufn_num = BUFSIZE-1;
  }
  return bufn++;
}


/* *****************************************************
   insert a suffix in the trie rooted at *p.
   we know that the trie already contains a string
   which share the first n chars with suf
   ***************************************************** */
static void insert_suffix(node *h, Int32 suf, int n, UChar mmchar)
{
   __inline__ node *new_node__blind_ssort(void);
  Int32 t;
  UChar c, *s;
  node *p, **pp;

  s = Text + suf;

#if 0
  // ---------- find the insertion point
  while( (t=h->skip) < n) {
    if( t < 0) break;          // insert "suf" just above *h
    c=s[t];  p=h->down;  
    // in the list p there must be a branch corresponding to c
  repeat:
    if(c==p->key) {              // found branch corresponding to c
      h = p;                     // go down
      continue;
    }
    if(c>p->key && p->right!=NULL) {
      p=p->right;
      goto repeat;
    }
    // no branch corresponding to c exists: this is a fatal error
    fprintf(stderr,"Error in blind_sort (insert_string)\n");
    exit(1);
  }
#else
  for(t=0;t<Stack_size;t++) {
    h=Stack[t];
    if(h->skip<0 || h->skip>=n) break;
  }  
#endif
  
  assert(s[n]!=mmchar || h->skip==-1 || h->skip==n);

  // --------- insert a new node before node *h if necessary
  if(h->skip!=n) {
    p = new_node__blind_ssort();     // create and init new node
    p->key = mmchar;
    p->skip = h->skip;  // p inherits skip and children of *h
    p->down = h->down;   
    p->right = NULL;
    h->skip = n;
    h->down = p;        // now *h has p as the only child 
  }
  assert(h->skip==n);

  // -------- search the position of s[n] among *h offsprings
  c=s[n]; pp = &(h->down);
  while((*pp)!=NULL) {
    if((*pp)->key>=c) 
      break;
    pp = &((*pp)->right);
  }
  // ------- insert new node containing suf
  p = new_node__blind_ssort();
  p->skip = -1;
  p->key = c; 
  p->right = *pp; *pp = p;
  p->down = (void *) suf;
  return;
}

/* ************************************************************
   this procedures traverse the trie in depth first order
   so that the suffixes (stored in the leaf) are recovered
   in lexicographic order
   ************************************************************ */
static void traverse_trie(node *h)
{
  node *p, *nextp;

  if(h->skip<0)
    Aux[Aux_written++] = (Int32) h->down;
  else {
    p = h->down;
    assert(p!=NULL);
    do {
      nextp = p->right;
      if(nextp!=NULL) {
	assert(nextp->key>=p->key);
	// if there are 2 nodes with equal keys 
	// they must be considered in inverted order
	if(nextp->key==p->key) {
	  traverse_trie(nextp);
	  traverse_trie(p);
	  p = nextp->right;
	  continue;
	}
      }
      traverse_trie(p);
      p=nextp;
    } while(p!=NULL);
  }
}



/* ***********************************************************************
   Function to compute the lcp of two strings originating from the *b1 and *b2
   the parameter is the length of s1 (which is shortest than s2)
   if s1 is a prefix of s2 we return the length of s1 -1
   The size of the unrolled loop must be at most equal to the costant 
   Cmp_overshoot defined in common.h
   the function return the result of the comparison (+ or -) and writes 
   in Cmp_done the number of comparisons done
   *********************************************************************** */ 
static __inline__
Int32 get_lcp_unrolled(UChar *b1, UChar *b2, Int32 cmp_limit)
{
  Int32 cmp2do; 
  UChar c1, c2;
  assert(b1 != b2);

  // execute blocks of 16 comparisons untill a difference
  // is found or we reach cmp_limit comparisons
  cmp2do = cmp_limit;
  do {
    // 1
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      break;}
    b1++; b2++; 
    // 2
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -=  1; break; }
    b1++; b2++; 
    // 3
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -=  2; break; }
    b1++; b2++; 
    // 4
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -=  3; break; }
    b1++; b2++; 
    // 5
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -=  4; break; }
    b1++; b2++; 
    // 6
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -=  5; break; }
    b1++; b2++; 
    // 7
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -=  6; break; }
    b1++; b2++; 
    // 8
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -=  7; break; }
    b1++; b2++; 
    // 9
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -=  8; break; }
    b1++; b2++; 
    // 10
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -=  9; break; }
    b1++; b2++; 
    // 11
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -= 10; break; }
    b1++; b2++; 
    // 12
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -= 11; break; }
    b1++; b2++; 
    // 13
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -= 12; break; }
    b1++; b2++; 
    // 14
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -= 13; break; }
    b1++; b2++; 
    // 15
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -= 14; break; }
    b1++; b2++; 
    // 16
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -= 15; break; }
    b1++; b2++; 

    cmp2do -= 16;
  } while(cmp2do>0);


  if(cmp_limit - cmp2do < cmp_limit)
    return cmp_limit-cmp2do;

  return cmp_limit-1;
} 



/* ************************************************************************
   this function returns the lcp between suf1 and suf2 (that is returns n 
   such that suf1[n]!=suf2[n] but suf1[i]==suf2[i] for i=0..n-1
   However, it is possible that suf1 is a prefix of suf2 (not vice-versa
   because of the initial sorting of suffixes in order of descreasing length)
   in this case the function returns n=length(suf1)-1. So in this case 
   suf1[n]==suf2[n] (and suf1[n+1] does not exists). 
   ************************************************************************ */
static Int32 compare_suffixes(Int32 suf1, Int32 suf2, Int32 depth)
{
  __inline__ Int32 get_lcp_unrolled(UChar *, UChar *, Int32);  
  int limit;
  UChar *s1, *s2;

  assert(suf1>suf2);
  s1  = Text + depth +suf1;
  s2  = Text + depth +suf2;
  limit = Text_size - suf1 - depth;
  return depth + get_lcp_unrolled(s1 ,s2, limit);
}


  
/* ****************************************************************** 
   comparison function used to sort suffixes in order of 
   increasing length. Since suffixes are represented by their offset
   in the array, we sort these offsets in order of decreasing length.
   ****************************************************************** */
static int neg_integer_cmp(const void *a, const void *b)
{
  return *((Int32 *) b) -  *((Int32 *) a); 
}


// free memory used for trie nodes
static void free_node_mem()
{
  int i;

  for(i=free_num-1;i>=0;i--) {
    assert(freearr[i]!=NULL);
    free(freearr[i]);
  }
  // clear counters
  bufn_num=free_num=0;
}
















