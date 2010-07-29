/**
 * @file hash_table_string.h
 *
 * This module implements the hash_table datatype for 
 * <key,value> pairs: <char*,STRING_LIST_T*> .
 *
 * $Id: hash_table_string.h 2739 2008-06-03 04:19:34Z f.buske $
*/

// hash_table_string.h 
//
// HASH TABLE object definitions
//
#ifndef HASH_H
#define HASH_H
#include "macros.h"
#ifdef DMALLOC
#include "dmalloc.h"
#endif
#include "string-list.h"

// Data structures

#define MAX_BINS 10000000 // 10 million

/* hash table entry */
typedef struct hash_table_entry HASH_TABLE_STRING_ENTRY; 
struct hash_table_entry {
  char *key1;				/* character key */ 
  int key2;					/* integer key */
  STRING_LIST_T* hash_value; // hash value of <keq1, key2>
  HASH_TABLE_STRING_ENTRY *next;	/* link to collision list */
  HASH_TABLE_STRING_ENTRY *prev;	/* backward link in collision list */
};

/* hash table */
typedef struct hash_table *HASH_TABLE_STRING;
struct hash_table {
  int n;					/* number of bins in hash table */
  int n_entries;            /* total number of entries in the hash table */
  HASH_TABLE_STRING_ENTRY **table;	/* array of entry pointers */
};

// External Functions
extern int get_num_entries(
  HASH_TABLE_STRING ht
);

HASH_TABLE_STRING hash_create(
  int n
);

void hash_destroy(
  HASH_TABLE_STRING ht			/* hash table to destroy */
);

void hash_entry_destroy(
  HASH_TABLE_STRING_ENTRY *hte         /* hash table entry to destroy */
);

BOOLEAN hash_insert_str(
  char *key1,			/* character key */
  HASH_TABLE_STRING ht			/* the hash table */
);

BOOLEAN hash_insert(
  char *key1,    	/* character key */
  int key2,			/* integer key */
  HASH_TABLE_STRING ht         /* the hash table */
);

BOOLEAN hash_remove_str(
  char *key1,			/* character key */
  HASH_TABLE_STRING ht			/* the hash table */
);

BOOLEAN hash_remove(
  char *key1,                   /* character key */
  int key2,                     /* integer key */
  HASH_TABLE_STRING ht                 /* the hash table */
);

HASH_TABLE_STRING_ENTRY *hash_lookup(
  char *key1,		/* character key */
  int key2, 		/* integer key */
  HASH_TABLE_STRING ht,         /* the hash table */
  STRING_LIST_T* hash_value	// hash value of keys
);

HASH_TABLE_STRING_ENTRY *hash_lookup_str(
  char *key1,			/* character key */
  HASH_TABLE_STRING ht,		/* the hash table */
  STRING_LIST_T* hash_value		// hash value of key
);

#endif 
