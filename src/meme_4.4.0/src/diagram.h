/*
 * $Id: diagram.h 1048 2006-07-06 20:07:44Z cegrant $
 * 
 * $Log$
 * Revision 1.1  2005/07/29 18:36:45  nadya
 * Initial revision
 *
 */

#ifndef DIAGRAM_H
#  define DIAGRAM_H

#include "macros.h"
#include "logodds.h"

DEXTERN(int, norder, 0);/* number of motifs for which order is given */
EXTERN int order[MAXLO];/* left-to-right order of motifs */
EXTERN int space[MAXLO];/* distance from end of previous motif to start this */

/* motif block diagram in mast format (see diagram.y):
        diagram         <- mdiagram | sdiagram
        mdiagram        <- motif | motif-diagram
        sdiagram        <- spacer | spacer-mdiagram
        motif           <- [integer]
        spacer          <- integer
*/
DEXTERN(char *, diagram, (char *)0);	
EXTERN int dptr;	/* current read position in diagram */

/* prevent compiler complaints */
#if defined(__STDC__) || defined(__cplusplus)
extern int yyparse(void);
#else
extern int yyparse();
#endif


#endif
