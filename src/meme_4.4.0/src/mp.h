/*
 * $Id: mp.h 1048 2006-07-06 20:07:44Z cegrant $
 * 
 * $Log$
 * Revision 1.1  2005/07/29 18:45:39  nadya
 * Initial revision
 *
 */

/***********************************************************************
*								       *
*	MEME++							       *
*	Copyright 1994, The Regents of the University of California    *
*	Author: Bill Grundy 					       *
*								       *
***********************************************************************/
#ifndef MPMYID_H
#define MPMYID_H

extern void mpInit(
  int *argc,
  char **argv[]
);
extern void mpFinalize();
extern int mpNodes();
extern int mpMyID();
extern void mpReduceAdd(int *data);
extern void mpBroadcast(void *data, int bytes, int broadcaster);

#endif
