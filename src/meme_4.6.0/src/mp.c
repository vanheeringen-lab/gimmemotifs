/*
 * $Id: mp.c 66 2005-07-29 17:15:48Z nadya $
 * 
 * $Log$
 * Revision 1.1  2005/07/29 17:18:06  nadya
 * Initial revision
 *
 */

/***********************************************************************
*								       *
*	MEME++							       *
*	Copyright 1994, The Regents of the University of California    *
*	Author: Timothy L. Bailey, Bill Grundy			       *
*								       *
***********************************************************************/
/* mp.c */
#ifdef PARALLEL

#include <mpi.h>
#include "stdio.h"

/***********************************************************************/
/*
	mpInit

	Wrapper on MPI_Init
*/
/***********************************************************************/
extern void mpInit(
  int *argc,
  char **argv[]
)
{
  MPI_Init(argc, argv);

  /*
    Associate an error handler with all communicators
    so that they return errors instead of aborting.
    MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
  */

#ifdef DEBUG
  { 
    int rank, namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(processor_name,&namelen);

#undef fprintf
    fprintf(stderr, "howdy %d %s\n", rank, processor_name);
  }
#endif
} /* mpInit */

/***********************************************************************/
/*
	MPI_Finalize

	Wrapper on MPI_Finalize
*/
/***********************************************************************/
extern void mpFinalize()
{
  MPI_Finalize();
}

/***********************************************************************
 * int mpNodes
 * 
 * Return the total number of processing nodes.
 ***********************************************************************/
int mpNodes()
{
  static int size;

  if (size == 0)
    MPI_Comm_size(MPI_COMM_WORLD, &size);
  return size;
}

/***********************************************************************
 * int mpMyID
 * 
 * Return the ID of the current node. Nodes are numbered from 0.
 ***********************************************************************/
int mpMyID()
{
  static int rank;

  if (rank == 0)
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
}

/***********************************************************************
 * void mpReduceAdd
 * 
 * Do an addition reduction on a single, integral value across all
 * processors and store the result at all processors.
 ***********************************************************************/
void mpReduceAdd(int *data)
{
  int send = *data;

  MPI_Allreduce(&send, data, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
}

/***********************************************************************
 * void mpBroadcast
 * 
 * Broadcast an untyped sequence of bytes to all nodes.
 ***********************************************************************/
void mpBroadcast(void *data, int bytes, int broadcaster)
{
  MPI_Bcast(data, bytes, MPI_BYTE, broadcaster, MPI_COMM_WORLD);
}

#endif /* PARALLEL */
