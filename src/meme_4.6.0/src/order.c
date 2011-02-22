/***********************************************************************
 * FILE: order.c
 * AUTHOR: William Stafford Noble
 * PROJECT: MHMM
 * COPYRIGHT: 2001-2008, WSN
 * DESCRIPTION: Data structure for representing the order and spacing
 *              of motifs in a linear model.
 ***********************************************************************/
#include <string.h>
#include <assert.h>
#include "order.h"
#include "string-list.h"

char* END_TRANSITION = "*";

/*
 * The motifs are indexed from 1 to n.  In motif_spacers, the value at
 * position i is the length of the spacer following motif i, so
 * position 0 contains the length of the spacer before the first
 * motif.
 * */
struct order {
  double         pvalue;         /* p-value of sequence from which this order
                  and spacing was derived. */
  int            num_occurs;     // Total number of motif occurrences. 
  int            num_alloc;      // Number of motif occurrence slots.
  int            num_distinct;   // Number of distinct motifs. 
  STRING_LIST_T* motif_occurs;   // Order of motif occurrences. 
  int*           motif_spacers;  // Distance between motif occurrences. 
};

/***********************************************************************
 * Does this order object have any space left in it?
 ***********************************************************************/
BOOLEAN_T is_full
  (ORDER_T* order)
{
  return(order->num_alloc == order->num_occurs);
}

/***********************************************************************
 * Create an order-and-spacing object from a user-specified string.
 ***********************************************************************/
ORDER_T* create_order
  (char* order_string)
{
  ORDER_T* order;     // The object being created. 
  int num_equals;     // Number of equals in the given string. 
  int i_string;       // Index into the given string. 
  int string_length;  // The length of the given string. 
  int i_spacer;       // Index of the current spacer. 
  int i_motif;        // Index of the current motif. 
  int start_position; // Index in the given string of start of this entry. 

  // Don't do anything if no string was given. 
  if (order_string == NULL) {
    return NULL;
  }

  // Allocate the data structure. 
  order = (ORDER_T*)mm_malloc(sizeof(ORDER_T));

  // Count the number of equals signs in the string. 
  string_length = strlen(order_string);
  num_equals = 0;
  for (i_string = 0; i_string < string_length; i_string++) {
    if (order_string[i_string] == '=') {
      num_equals++;
    }
  }
  order->num_alloc = (int)(num_equals / 2);
  order->num_occurs = 0;

  // Allocate the arrays.
  order->motif_occurs = new_string_list();
  order->motif_spacers = (int*)mm_calloc((order->num_alloc)+2, sizeof(int));

  // Step through the string, looking for equals signs. 
  start_position = 0;
  i_spacer = 0;
  i_motif = 1;
  for (i_string = 0; i_string < string_length; i_string++) {
    if (order_string[i_string] == '=') {
      order_string[i_string] = '\0';
      if (i_spacer < i_motif) {
        if (sscanf(&(order_string[start_position]), "%d", 
            &(order->motif_spacers[i_spacer])) != 1) {
              die("Error reading order and spacing at spacer %d.\n", i_spacer);
  }
  if (order->motif_spacers[i_spacer] == 0) {
    die("Sorry, Meta-MEME cannot represent spacers of length 0.");
  }
  i_spacer++;
      } else {
  char motif_id[MAX_MOTIF_ID_LENGTH + 1];

  if (sscanf(&(order_string[start_position]), "%s", motif_id) != 1) {
    die("Error reading order and spacing at motif %d.\n", i_motif);
  }

  // Find out if this is a new motif.
  if (!order_contains(motif_id, order)) {
    order->num_distinct++;
  }

  // Add the motif ID to the list.
  add_string(motif_id, order->motif_occurs);
  order->num_occurs++;
  i_motif++;
      }
      start_position = i_string + 1;
    }
  }
  if (sscanf(&(order_string[start_position]), "%d", 
       &(order->motif_spacers[i_spacer])) != 1) {
    die("Error reading order and spacing at spacer %d.\n", i_spacer);
  }
  assert(is_full(order));

  // Tell the user. 
  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "Using order: ");
    print_order_and_spacing(stderr, order);
  }

  return(order);
}

/***********************************************************************
 * Create an empty order-and-spacing object.
 ***********************************************************************/
ORDER_T* create_empty_order
  (int    num_alloc,
   double sequence_p)
{
  ORDER_T* new_order;

  // Allocate the data structure. 
  new_order = (ORDER_T*)mm_malloc(sizeof(ORDER_T));

  new_order->num_alloc = num_alloc;
  new_order->num_occurs = 0;
  new_order->pvalue = sequence_p;
  new_order->num_distinct = 0;

  // Allocate the arrays. 
  new_order->motif_occurs = new_string_list();
  // Add 1 for end transition.
  new_order->motif_spacers = (int*)mm_calloc(num_alloc + 1, sizeof(int));

  return new_order;
}

/***********************************************************************
 * Add one spacer and the following motif to an order object.
 *
 * Does nothing if the order parameter is NULL.
 ***********************************************************************/
void add_occurrence
  (char*    motif_id,
   int      spacer_length,
   ORDER_T* order)
{
  extern char* END_TRANSITION;

  if (order == NULL) {
    return;
  }

  // Count distinct motifs.
  if (!order_contains(motif_id, order)) {
    order->num_distinct++;
  }

  // Record the spacer before this occurrence.
  order->motif_spacers[order->num_occurs] = spacer_length;

  // Don't include the end transition as a motif occurence
  if (strncmp(motif_id, END_TRANSITION, MAX_MOTIF_ID_LENGTH) != 0) {
    (order->num_occurs)++;
    add_string(motif_id, order->motif_occurs);
  }
}

/***********************************************************************
 * Determine whether a given motif order-and-spacing object contains a
 * given motif ID.
 ***********************************************************************/
BOOLEAN_T order_contains
  (char*    motif_id,
   ORDER_T* order)
{
  int i_order;

  for (i_order = 0; i_order < order->num_occurs; i_order++) {
    if (strcmp(get_nth_string(i_order, order->motif_occurs), motif_id) == 0) {
      return(TRUE);
    }
  }
  return(FALSE);
}

/***********************************************************************
 * Get the length of the spacer after a given motif.  Motifs are
 * indexed from 1, so 0 gets the initial spacer.
 ***********************************************************************/
int get_spacer_length
  (ORDER_T* order,
   int      motif_index)
{
  return(order->motif_spacers[motif_index]);
}

/***********************************************************************
 * Get and set the p-value of a given order object.
 *
 * The p-value of a null object is 1.0.
 ***********************************************************************/
void set_pvalue
  (double pvalue,
   ORDER_T* order)
{
  assert(order != NULL);
  order->pvalue = pvalue;
}

double get_pvalue
  (ORDER_T* order)
{
  if (order == NULL) {
    return(1.0);
  }
  return order->pvalue;
}

/***********************************************************************
 * Get the number of distinct motif occurrences.
 ***********************************************************************/
int get_num_distinct
  (ORDER_T* order)
{
  if (order == NULL) {
    return(0);
  }
  return order->num_distinct;
}


/***********************************************************************
 * Re-order the motifs in the order in which they will appear in the model.
 ***********************************************************************/
void reorder_motifs
  (ORDER_T*  order,
   int*      num_motifs,
   MOTIF_T*  motifs)
{
  MOTIF_T* temp_motifs; /* Pointer to a copying buffer. */
  int      i_motif;     /* Index of the desired motif. */
  int      i_order;     /* Index into the order object. */

  assert(order != NULL);
  
  /* Allocate memory for a temporary copy of the motif array. */
  temp_motifs = (MOTIF_T *)mm_malloc(sizeof(MOTIF_T) * order->num_alloc);

  /* Print the original list. */
  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "Motifs:");
    for (i_motif = 0; i_motif < *num_motifs; i_motif++) {
      fprintf(stderr, " %s", motifs[i_motif].id);
    }
    fprintf(stderr, ".\n");
    fprintf(stderr, "Motif order: ");
    print_order_and_spacing(stderr, order);
  }

  /* Copy the original motifs into the temporary space in the right order. */
  for (i_order = 0; i_order < order->num_occurs; i_order++) {
    BOOLEAN_T found_motif_token = FALSE;
    char *order_motif_token = get_nth_string(i_order, order->motif_occurs);
    for (i_motif = 0; i_motif < *num_motifs; i_motif++) {
      if (strcmp(motifs[i_motif].id, order_motif_token) == 0) {
        found_motif_token = TRUE;
        copy_motif(&(motifs[i_motif]), &(temp_motifs[i_order]));
        break;
      }
    }
    if (found_motif_token == FALSE) {
      die(
          "The motif %s in the order string was not "
          "found in the set of allowed motifs.",
          order_motif_token
      );
    }
  }
  // Free the original motifs since we've copied them 
  // and are about to write over them. (free_motif actually frees
  // the frequency matrix contained in the motif).
  for (i_motif = 0; i_motif < *num_motifs; i_motif++) {
    free_motif(&(motifs[i_motif])); 
  }

  /* Copy back into the original place. */
  for (i_order = 0; i_order < order->num_occurs; i_order++) {
    copy_motif(&(temp_motifs[i_order]), &(motifs[i_order]));
    free_motif(&(temp_motifs[i_order]));
  }
  *num_motifs = order->num_occurs;

  /* Print the re-ordered list. */
  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "Re-ordered motifs: ");
    for (i_motif = 0; i_motif < *num_motifs; i_motif++) {
      fprintf(stderr, "%s ", motifs[i_motif].id);
    }
    fprintf(stderr, "\n");
  }

  /* Free the temporary space. */
  myfree(temp_motifs);
}

/***********************************************************************
 * Print the order and spacing with equal signs between.
 ***********************************************************************/
void print_order_and_spacing
  (FILE*    outfile,
   ORDER_T* order)
{
  int i_order;

  fprintf(outfile, "%d", order->motif_spacers[0]);
  for (i_order = 0; i_order < order->num_occurs; i_order++) {
    fprintf(outfile, "=%s=%d", get_nth_string(i_order, order->motif_occurs),
      order->motif_spacers[i_order + 1]);
  }
  fprintf(outfile, "\n");
}

/***********************************************************************
 * Free dynamic memory used by an order object.
 ***********************************************************************/
void free_order
  (ORDER_T* order)
{
  /* Don't deallocate an empty struct. */
  if (order == NULL) {
    return;
  }

  free_string_list(order->motif_occurs);
  myfree(order->motif_spacers);
  myfree(order);
}
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

