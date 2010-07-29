/*
 * @file seed.h 
 *
 * A partition object, used for indicating blocks within the matrix of an
 * SP_MATRIX object.
 *
 */

#ifndef PARTITION_H
#define PARTITION_H

#include "macros.h"
#ifdef DMALLOC
#include "dmalloc.h"
#endif


/**
 * A PARTITION type.
 */


/**
 * A PARTITION type.
 */
typedef struct partition PARTITION;
struct partition {
  int min_w;
  int max_w;
  int central_w;
  int min_n;
  int max_n;
  int central_n;
};


extern PARTITION *new_partition(
  int part_min_w,   ///< Minimum w value in the partition
  int part_max_w,   ///< Maximum w value in the partition
  int central_w,    ///< Central w value of this partition
  int part_min_n,   ///< Minimum n value in the partition
  int part_max_n,   ///< Maximum n value in the partition
  int central_n     ///< Central w value of this partition
);


extern void print_partition (
  PARTITION *part,   ///< The partition being printed
  FILE *out          ///< The output destination.
);

#endif


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
