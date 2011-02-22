#####################################################################
## FILE: Makefile.am
# AUTHOR: William Stafford Noble, Charles E. Grant
# CREATE DATE: 2003
# PROJECT: Meta-MEME
# COPYRIGHT: 2003, UW
#####################################################################

## Automake uses this file to generate a Makefile

# This Makefile builds the test output used to test FIMO.

# These files should be removed as part of the 'make clean' target.
CLEANFILES = 

# These files should be included with the distribution.
DATAFILES = \
  fimo.txt \
  fimo-motif23.txt \
  fimo-priors.txt \
  GCN4.meme \
  GCN4_YPD.fasta \
  GCN4_YPD.psp \
  MCM1.meme.html \
  prior.dist.txt \
  spiked.fasta

# If .svnignore file doesn't exist
# create an empty one
.svnignore:
	if [ -e .svn ] ; then touch .svnignore ; fi

# Update SVN svn:ignore property
# If .svnignore exists and longer then 0 bytes.
.svnignore.old:.svnignore
	if [ -e .svn ] ; \
  then \
    if [ -s .svnignore ] ; \
    then \
      svn propset svn:ignore -F .svnignore .  ; \
      cp .svnignore .svnignore.old; \
    fi \
  fi


EXTRA_DIST = $(DATAFILES)
