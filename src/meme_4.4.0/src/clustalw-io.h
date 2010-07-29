/****************************************************************************
 * FILE: clustalw-io.h
 * AUTHOR: William Stafford Noble, Charles E. Grant
 * CREATE DATE: 10/27/2004
 * PROJECT: EVOMCAST
 * DESCRIPTION: Read an alignment from a CLUSTALW file (.aln) into memory.
 * COPYRIGHT: 1998-2008 WSN
 ****************************************************************************/
#ifndef CLUSTALW_IO_H
#define CLUSTALW_IO_H

#include "alignment.h"

#define BLOCKSIZE 60

/****************************************************************************
 * Read an alignment from an open FILE in CLUSTALW (.aln) format.
 * If call is successful caller is responsible for freeing alignment object.
 *
 * Return: Was an alignment successfully read?
 ****************************************************************************/
BOOLEAN_T read_clustalw (FILE* clustalw_file, ALIGNMENT_T** alignment);

/****************************************************************************
 * Create an alignment obbject from the named file file in CLUSTALW 
 * (.aln) format. If call is successful, caller is responsible for 
 * freeing alignment object.
 *
 * Return: The alignment that was read from the file.
 ****************************************************************************/
ALIGNMENT_T* read_alignment_from_clustalw_file(char* filename);

/****************************************************************************
 * Print an alignment in CLUSTALW format.
 ****************************************************************************/
void print_clustalw( 
    FILE* outfile,
    BOOLEAN_T show_residue_count,
    ALIGNMENT_T* alignment
);

#endif
