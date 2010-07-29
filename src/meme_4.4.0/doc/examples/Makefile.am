## $Id: Makefile.am 4514 2010-05-04 02:28:44Z james_johnson $
##
## Process this file with automake to produce Makefile.in
##
## # $Log$
## # Revision 1.3  2005/09/15 21:35:28  nadya
## # remove extra trailing slash
## #
## # Revision 1.2  2005/09/14 00:39:52  nadya
## # add mast sample sequences
## #
## # Revision 1.1.1.1  2005/07/25 23:21:55  nadya
## # Importing from meme-3.0.14, and adding configure/make
## #
##

AUTOMAKE_OPTIONS = no-dependencies

SUBDIRS = \
  fimo_example_output_files \
  glam2_example_output_files \
  glam2scan_example_output_files \
  meme_example_output_files \
	mast_example_output_files \
	mcast_example_output_files \
	tomtom_example_output_files

EXAMPLES = \
	At.fa \
  sample-kabat.seq \
  sample-glam2-aln.dna \
  sample-glam2-aln.protein \
  sample-dna-motif.meme-io \
  sample-dna-motif.mast-io \
  sample-protein-motif.meme-io \
  sample-protein-motif.mast-io

if WEBSITE
exdir = $(WEBDIR)/examples
else
exdir = $(prefix)/doc/examples
endif

ex_DATA = $(EXAMPLES)


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

EXTRA_DIST = $(EXAMPLES)
