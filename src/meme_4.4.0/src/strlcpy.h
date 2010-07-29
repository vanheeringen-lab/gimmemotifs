#ifndef STRLCPY_H_
# define STRLCPY_H_ 1

# if HAVE_CONFIG_H
#  include <config.h>
# endif

# ifndef PARAMS
#  if defined PROTOTYPES || (defined __STDC__ && __STDC__)
#   define PARAMS(Args) Args
#  else
#   define PARAMS(Args) ()
#  endif
# endif

#include <stdio.h>

size_t strlcpy PARAMS ((char *, const char *, size_t));

#endif  /* STRLCPY_H_ */

