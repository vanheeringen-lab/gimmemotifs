## $Id: Makefile.am 1421 2006-11-01 01:15:40Z cegrant $
##
## Process this file with automake to produce Makefile.in

AUTOMAKE_OPTIONS = no-dependencies

DATAFILES = \
	all-A.fasta \
	one-peak-dna.fasta \
	one-peak-protein.fasta \
	README \
	one-peak-dna-revcomp.psp \
	one-peak-protein.psp

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
