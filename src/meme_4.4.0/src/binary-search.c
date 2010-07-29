#include "binary-search.h"

#define GET(INDEX) ((void*)(((char*)base)+(size*(INDEX))))
/*
 * binary_search
 * returns the index of the item equal to the key as compared using the compar function. If an equal item does not exist
 * then return -(indx + 1) where indx is the index that the item should be inserted. If the list is not sorted then the result
 * is undefined though it will return something...
 * The compar function works exactly like in bsearch so look that up if you need more details.
 */
int binary_search(const void *key, const void *base, size_t nmemb, size_t size, int (*compar)(const void *, const void *)) {
  int low, high, mid, cmp;
  low = 0;
  high = nmemb;
  while (low < high) {
    mid = low + ((high - low) / 2);
    if (compar(key, GET(mid)) > 0) {
      low = mid + 1;
    } else {
      high = mid;
    }
  }
  if ((low < nmemb) && (compar(key, GET(low)) == 0)) {
    return low;
  } else {
    return -(low + 1);
  }
}
