/*
 * $Id: aa.h 1048 2006-07-06 20:07:44Z cegrant $
 * 
 * $Log$
 * Revision 1.1  2005/07/29 18:34:01  nadya
 * Initial revision
 *
 */

/**********************************************************************/
/*
	associative array data structures
*/
/**********************************************************************/

#ifndef AA_H
#define AA_H

#include <macros.h>

#define AA_BUFSIZE 1000

/* types */

/* names of associative array value types */
typedef enum {Int, Double, String} TYPE_NAME;

/* union of types for values in associative arrays */
typedef union value VALUE;
union value {
 int i; 
 double d;
 char *s;
};

/* associative array entry */
typedef struct aa_entry AA_ENTRY; 
struct aa_entry {
  char *key;			/* key */ 
  TYPE_NAME value_type;		/* type of value */
  VALUE value;			/* value */
  AA_ENTRY *next;		/* link to collision list */
};

/* associative array */
typedef struct aa *AA;
struct aa {
  int n;			/* number of bins in hash table */
  AA_ENTRY **table;		/* array of entry pointers */
};


/* functions */

extern AA aa_create(
  int n					/* size of associative array */
);

extern void aa_destroy(
  AA aa 				/* associative array to destroy */
);

extern AA_ENTRY *aa_lookup(
  AA aa,                                /* the associative array */
  char *key                             /* key */
);

extern void aa_assign(
  AA aa,				/* the associative array */
  char *key,				/* key */
  TYPE_NAME value_type,			/* type of value to assign */
  VALUE value				/* value to assign */
);

extern char **aa_keys(
  AA aa					/* the associative array */
);

/* macros */

/* look up the value of a given type in an associative array */
#define aa_get_int(aa, key) ( (aa_lookup((aa), (key)))->value.i )
#define aa_get_double(aa, key) ( (aa_lookup((aa), (key)))->value.d )
#define aa_get_string(aa, key) ( (aa_lookup((aa), (key)))->value.s )

/* assign the value of a given type to an associative array key */
#define aa_put_int(aa, key, value) {					\
  VALUE val; 								\
  val.i = value;							\
  aa_assign(aa, key, Int, val);						\
}
#define aa_put_double(aa, key, value) {					\
  VALUE val; 								\
  val.d = value;							\
  aa_assign(aa, key, Double, val);					\
}
#define aa_put_string(aa, key, value) {					\
  VALUE val; 								\
  val.s = value;							\
  aa_assign(aa, key, String, val);					\
}

/* see if a key is defined in an associative array */
#define aa_defined(aa, key) ( aa_lookup((aa), (key)) ? TRUE : FALSE )

#endif

