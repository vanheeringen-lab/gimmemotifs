if BUILD_LIBXSLT
noinst_LTLIBRARIES = libxslt.la

libxslt_la_SOURCES = 			\
	attrvt.c			\
	xslt.c				\
	xsltutils.c			\
	pattern.c			\
	templates.c			\
	variables.c			\
	keys.c				\
	numbers.c			\
	extensions.c			\
	extra.c				\
	functions.c			\
	namespaces.c			\
	imports.c			\
	attributes.c			\
	documents.c			\
	preproc.c			\
	transform.c			\
	security.c			\
	win32config.h			\
	xsltwin32config.h		\
	xsltwin32config.h.in		\
	xslt.h				\
	xsltutils.h			\
	pattern.h			\
	templates.h			\
	variables.h			\
	keys.h				\
	numbersInternals.h		\
	extensions.h			\
	extra.h				\
	functions.h			\
	namespaces.h			\
	imports.h			\
	attributes.h			\
	documents.h			\
	preproc.h			\
	transform.h			\
	security.h			\
	xsltInternals.h			\
	xsltconfig.h			\
	xsltexports.h \
	libxslt.h

libxslt_la_CFLAGS = $(AM_CFLAGS) $(LIBXML2_CFLAGS) -I$(top_srcdir)/src
libxslt_la_LIBADD = $(LIBXML2_LIBS)

EXTRA_DIST = Copyright

cleanup:
	-@(find . -name .\#\* -exec rm {} \;)

endif

