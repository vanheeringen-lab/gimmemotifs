/*
 * @file partition.c 
 *
 * A partition object, used for indicating blocks within the matrix of an
 * SP_MATRIX object.
 * 
 */

#include "macros.h"
#include "partition.h"
#ifdef DMALLOC
#include "dmalloc.h"
#endif


/**
 * new_partition generates a new PARTITION object given the specified
 * (inclusive) borders of the partition
 *
 * \return A pointer to the new object.
 */
extern PARTITION *new_partition(
  int part_min_w,   ///< Minimum w value in the partition
  int part_max_w,   ///< Maximum w value in the partition
  int central_w,    ///< Central w value of this partition
  int part_min_n,   ///< Minimum n value in the partition
  int part_max_n,   ///< Maximum n value in the partition
  int central_n     ///< Central w value of this partition
)
{
  PARTITION *new_part = NULL;
  //new_part = (PARTITION *)mymalloc(sizeof PARTITION);
  Resize(new_part, 1, PARTITION);

  new_part->min_w = part_min_w;
  new_part->max_w = part_max_w;
  new_part->min_n = part_min_n;
  new_part->max_n = part_max_n;
  new_part->central_w = central_w;
  new_part->central_n = central_n;

  return(new_part);
}


/**
 * print_partition
 *
 * Print a representation of this partition.
 *
 */
extern void print_partition (
  PARTITION *part,   ///< The partition being printed
  FILE *out          ///< The output destination.
) {
  fprintf(out, "-------------PARTITION OBJECT-------------\n");
  fprintf(out, "min_w = %d, max_w = %d. min_n = %d, max_n = %d.\n",
          part->min_w, part->max_w, part->min_n, part->max_n);
  fprintf(out, "-----------END PARTITION OBJECT-----------\n");
} // print_partition


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
