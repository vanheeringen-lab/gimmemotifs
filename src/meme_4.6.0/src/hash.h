/*
 * $Id: hash.h 776 2006-05-10 17:27:13Z cegrant $
 * 
 * $Log$
 * Revision 1.1  2005/07/29 18:39:19  nadya
 * Initial revision
 *
 */

/**********************************************************************/
/*
	hash table data structures
*/
/**********************************************************************/
/* 5-24-00 tlb; fix hashing function and make generic */

#ifndef HASH_H
#define HASH_H

#include <macros.h>

/* hash table entry */
typedef struct hash_table_entry HASH_TABLE_ENTRY; 
struct hash_table_entry {
  char *key1;			/* character key */ 
  int key2;			/* integer key */
  HASH_TABLE_ENTRY *next;	/* link to collision list */
};

/* hash table */
typedef struct hash_table *HASH_TABLE;
struct hash_table {
  int n;			/* number of bins in hash table */
  HASH_TABLE_ENTRY **table;	/* array of entry pointers */
};

extern int hash(
  char *key1, 			/* character key */
  int key2,			/* integer key */
  int n				/* modulo */
);

extern HASH_TABLE hash_create(
  int n
);

extern void hash_destroy(
  HASH_TABLE ht			/* hash table to destroy */
);

extern void hash_insert(
  char *key1,    	/* character key */
  int key2,		/* integer key */
  HASH_TABLE ht         /* the hash table */
);

extern BOOLEAN hash_lookup(
  char *key1,		/* character key */
  int key2, 		/* integer key */
  HASH_TABLE ht         /* the hash table */
);

#endif 

