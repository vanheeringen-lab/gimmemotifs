#####################################################################
## FILE: Makefile.am
# AUTHOR: William Stafford Noble, Charles E. Grant
# CREATE DATE: 2003
# PROJECT: Meta-MEME
# COPYRIGHT: 2003, UW
#####################################################################

## Automake uses this file to generate a Makefile

# This Makefile builds the scripts used to test the Meta-MEME components.

.PHONY: profile
profile:
	pushd ../../src; \
  $(MAKE) clean; \
  $(MAKE) AM_CFLAGS='$(AM_CFLAGS) -ggdb -O0 -pg  '; \
  popd; \
  ./cmd;

# These files should be removed as part of the 'make clean' target.
CLEANFILES = 

# These files should be included with the distribution.
DATAFILES = \
  MCM1.meme.html \
  motiph.gff \
  motiph-motif23.gff \
  shadow.gff \
  spiked.aln \
  spiked.fasta \
  yeast.tree

EXTRA_DIST = $(DATAFILES)
