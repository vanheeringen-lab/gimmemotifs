/********************************************************************
 * FILE: gomo_highlight.h
 * AUTHOR: James Johnson 
 * CREATE DATE: 02/10/2009
 * PROJECT: MEME suite
 * COPYRIGHT: 2009, UQ
 *
 * GOMO highlight uses the GO Hierarchy to mark results
 * as being implied when a more specific term in the hierarchy is
 * also a significant result. It also includes the grouping and name 
 * of the term and the relative position in the hierarchy of the term.
 *
 ********************************************************************/


#ifndef GOMO_HIGHLIGHT_H
#define GOMO_HIGHLIGHT_H


typedef struct GODAG GODAG_T;

/*
 * loads a gotree from file
 */
GODAG_T* load_go_dag(char *file);

/*
 * destroys a gotree
 */
void destroy_go_dag(GODAG_T *godag);

/*
 * Rewrites the gomo xml to include extra detail concerning the go hierarchy 
 */
void rewrite_gomo_xml(const char *tmp_dir, char *xml_file, GODAG_T *godag); 

#endif
