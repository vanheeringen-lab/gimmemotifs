## $Id: Makefile.am 4334 2010-01-18 09:13:21Z cegrant $
##
## Process this file with automake to produce Makefile.in
##

#
# To MAKE in DEBUG mode (it will write out lots of status info):
#	make CFLAGS=" -DDEBUG_PARALLEL"

# redefine make macros 
#CC = $(MPICC)
# TLB: FIXME This command makes things at least compile with the -DPARALLEL switch,
# but there are two -O switches
#CFLAGS=
#INCLUDES = -I.. $(MPIINC)
#COMPILE = $(CC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES)
#LINK = $(CCLD) $(AM_CFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@
#LINK = $(LIBTOOL) --mode=link --tag=CC $(CCLD) $(AM_CFLAGS) $(CFLAGS) \
#	$(AM_LDFLAGS) $(LDFLAGS) -o $@ 

CC = $(MPICC)

bin_PROGRAMS = meme_p
meme_p_LDADD = ../libcommon.la $(LIBXML2_LIBS) $(LIBXSLT_LIBS)
meme_p_CFLAGS = -I.. $(MPIFLAGS) $(LIBXML2_CFLAGS) $(LIBXSLT_CFLAGS)

meme_p_SOURCES = \
  ../banner.c \
  ../branching_search.c \
  ../calculate_p_y.c \
  ../discretize.c \
  ../display.c \
  ../dpalign.c \
  ../em.c  \
  ../init.c  \
  ../justlike.c \
  ../likelihood.c \
  ../llr.c \
  ../logs.c \
  ../meme.c \
  ../meme-print-html.c \
  ../meme_util.c \
  ../message.c \
  ../mp.c \
  ../oops.c \
	../partition.c \
  ../prior.c  \
  ../psp.c \
  ../pssm-distr.c  \
  ../read_seq_file.c \
  ../regress.c \
  ../seed.c \
  ../seed_diffs.c \
  ../seq2theta.c \
  ../sp_matrix.c \
  ../starts.c \
  ../subseq7.c \
  ../tcm.c \
  ../verbosity.c

