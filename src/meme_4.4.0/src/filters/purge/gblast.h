/* gblast.h - generic blast program. */
#if !defined (GBLAST)
#define GBLAST
#include <time.h>
#include "afnio.h"
#include "random.h"
#include "residues.h"
#include "alphabet.h"
#include "dheap.h"
#include "sequence.h"
#include "mlist.h"

/*********************** Finite State Machine ************************
FSM = (Q,q0,A,Sigma,d)
	Q = States
	q0 e Q = start state
	A subset Q = accepting states
	Sigma = input alphabet
	d = function from Q x Sigma -> Q (transition function)

if in an accepting state then execute appropriate action.
	- go to positions on list and try extending alignment

input text T[1..n];  pattern to be found P[1..m]
 n = input string length; m = pattern length 

	(note: lex is a FSM)

  input tokens = A-Y and '-' 'x' '*' '\0'

 if q = a then  go to query sequence:
   pos[q][1..n][A...Y] = list of positions matching pattern in accepting 
	state = NULL if not an accepting state.

  blast method:
	1. compile list of high scoring words and make fsm.
	2. scan database for hits.
	3. extend hits.
(for purge extend only until find that score >= cutoff.)

	QWL
	..:  -> S = R(Q,N) + R(W,Y) + R(L,L).
	NYL	
		if(S > Threshold) then extend hit to find MSP.

		need drop-score.
 *********************************************************************/

/*************************** generic gblast type **************************/
typedef struct {
	long	nQ;		/** number of States **/
	long	**d;		/** d[q][r] = function from Q x A -> Q **/
	ml_type	*pos;		/* lists for accept */
	long	*tmp;		/* temporary list */
	long	T;
	a_type  A;		/* alphabet */
	e_type  E;              /* query sequence */
} gblast_type;
typedef gblast_type *gb_typ;
/*********************************************************************/
/******************************* private *******************************/
Boolean FastExtendGBlast(e_type E1, long i1, e_type E2, long i2,
        register int **R, long score);
long     gblast_error(char *s);

/******************************* PUBLIC *******************************/
gb_typ	MakeGBlast(long T, e_type E, a_type A) ;
long     MatcherGBlastOffset(e_type E, gb_typ B, long *os);
void    NilGBlast(gb_typ B);
long     MatcherGBlast(FILE *fptr, e_type E, gb_typ B);
Boolean FastMatcherGBlast(e_type E, gb_typ B, long score);
long     ExtendGBlast(e_type E1, long i1, e_type E2, long i2, a_type A);
/*********************************************************************/

/* CODES */

/* CONSTANTS */
#define	USAGE2_START	"USAGE: fsm seq1 seq2 [options]\n\
   options:\n\
     [-l]	      - eliminate low complexity sequences\n\
     [-f]	      - create output file for sequence with repeats\n\
     [-F]	      - create output file for sequence without repeats\n\
     [-e<float>]      - E value cutoff for output of sequences\n"

#endif

