/*
 * $Id: hashu.h 1048 2006-07-06 20:07:44Z cegrant $
 * 
 * $Log$
 * Revision 1.1  2005/07/29 18:39:47  nadya
 * Initial revision
 *
 */

/**********************************************************************/
/*
	hash table data structures
*/
/**********************************************************************/

#ifndef HASHU_H
#define HASHU_H

#include <macros.h>

/* hash table value */
typedef struct hashu_value HASHU_VALUE;
struct hashu_value {
  int type;			/* type of entry: 0=char*, 1=int, 2=double */
  union {			/* actual value depending on type */
    char *cval;
    int ival;
    double dval;
  } data;
};

/* hash table entry */
typedef struct hashu_table_entry HASHU_TABLE_ENTRY; 
struct hashu_table_entry {
  char *key;			/* character key */ 
  HASHU_VALUE value;		/* value of entry */
  HASHU_TABLE_ENTRY *next;	/* link to collision list */
};

/* hash table */
typedef struct hashu_table *HASHU_TABLE;
struct hashu_table {
  int n;			/* number of bins in hash table */
  BOOLEAN nstrings;		/* number of strings in table */
  HASHU_TABLE_ENTRY **table;	/* array of entry pointers */
};

extern int hashu(
  char *key, 				/* key */
  int w,				/* key width; if 0, null-terminated */
  int n					/* size of hash table */
);

extern HASHU_TABLE hashu_create(
  int n					/* size of hash table */
);

extern void hashu_destroy(
  HASHU_TABLE ht 		/* hash table to destroy */
);

extern BOOLEAN hashu_insert(
  char *key, 			/* key */
  int w,			/* key width; if 0, null-terminated */
  HASHU_VALUE value,		/* value to insert */
  HASHU_TABLE ht		/* the hash table */
);

extern HASHU_VALUE hashu_lookup(
  char *key,			/* key */
  int w,			/* key width; if 0, null-terminated */
  HASHU_TABLE ht		/* the hash table */
);

extern double hashu_increment(
  char *key, 			/* key */
  int w,			/* key width; if 0, null-terminated */
  HASHU_TABLE ht		/* the hash table */
);

#endif

/* $$ */
