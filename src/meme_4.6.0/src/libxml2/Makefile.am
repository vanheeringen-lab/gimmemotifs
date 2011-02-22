## Process this file with automake to produce Makefile.in

if BUILD_LIBXML2
SUBDIRS = include .

DIST_SUBDIRS = include .
noinst_LTLIBRARIES = libxml2.la

INCLUDES = -I@srcdir@ -I@srcdir@/include -I@srcdir@/include/libxml2

# We don't control this code so turn off warnings.
libxml2_la_CFLAGS =  $(AM_CFLAGS) -w
libxml2_la_LIBADD =  $(ICONV_LIBS) 

libxml2_la_SOURCES = entities.c encoding.c error.c parserInternals.c  \
		parser.c tree.c hash.c list.c xmlIO.c xmlmemory.c uri.c  \
		valid.c HTMLparser.c HTMLtree.c debugXML.c xpath.c  \
		globals.c threads.c xmlstring.c \
		xmlregexp.c xmlschemas.c xmlschemastypes.c xmlunicode.c \
		xmlreader.c relaxng.c dict.c SAX2.c \
		chvalid.c pattern.c xmlsave.c \
		trionan.c

EXTRA_DIST = \
  Copyright \
  elfgcchack.h \
  libxml.h \
  triodef.h \
  trio.h \
  trionan.h \
  triop.h \
  triostr.h

DEPS = @srcdir@/libxml2.la

cleanup:
	-@(find . -name .\#\* -exec rm {} \;)
endif
