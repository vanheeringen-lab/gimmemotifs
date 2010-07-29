AC_PREREQ(2.59)
AC_INIT([gadem], [1.0], [li3@niehs.nih.gov])
AM_INIT_AUTOMAKE([-Wall foreign])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
Makefile
src/Makefile
examples/Makefile
])


# Checks for programs.
AC_PROG_CC
CFLAGS=""

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [pow])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([limits.h malloc.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T

AC_CHECK_SIZEOF(long)

AM_CONDITIONAL(LARGERINT, [test "$ac_cv_sizeof_long" = "8"])

# Checks for library functions.
AC_CHECK_FUNCS([pow fabs qsort floor rint])
AC_FUNC_MALLOC
AC_FUNC_REALLOC

AC_OUTPUT
