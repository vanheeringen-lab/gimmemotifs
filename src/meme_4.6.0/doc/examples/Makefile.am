## $Id: Makefile.am 5254 2010-12-23 05:35:13Z phillip_machanick $
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

##############################################
# some examples are also used in tests: check
# scripts/metatest.pl.in before changing

AUTOMAKE_OPTIONS = no-dependencies


SUBDIRS = \
  compute_prior_dist_example_output_files \
  fimo_example_output_files \
  glam2_example_output_files \
  glam2scan_example_output_files \
  meme_example_output_files \
  mast_example_output_files \
  mcast_example_output_files \
  sample_opal_scripts \
  tomtom_example_output_files

EXAMPLES = \
  At.fa \
  Jaspar-subset.meme \
  sample-kabat.seq \
  sample-glam2-aln.dna \
  sample-glam2-aln.protein \
  sample-dna-Klf1-200-100.fa \
  sample-dna-Klf1-200-100-ame.html \
  sample-dna-Klf1-200-100-dreme.txt \
  sample-dna-Klf1-200-100-ame.txt \
  sample-dna-Klf1-200-100-shuffled.fa \
  sample-dna-Klf1-200.fa \
  sample-dna-Klf1.fa \
  sample-dna-motif.meme-io \
  sample-dna-motif.mast-io \
  sample-protein-motif.meme-io \
  sample-protein-motif.mast-io

EXAMPLES_TAR = \
  memechip_example_output_files.tgz

if WEBSITE
exdir = $(WEBDIR)/examples
other_exdir = $(WEBDIR)/doc/examples
else
exdir = $(prefix)/doc/examples
other_exdir = $(exdir)
endif

ex_DATA = $(EXAMPLES)

install-data-hook:
	tar xzvf memechip_example_output_files.tgz -C $(other_exdir)

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

EXTRA_DIST = $(EXAMPLES) $(EXAMPLES_TAR)
