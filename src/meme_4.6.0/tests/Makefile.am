## $Id: Makefile.am 5207 2010-12-16 06:33:29Z phillip_machanick $
##
## Process this file with automake to produce Makefile.in

AUTOMAKE_OPTIONS = no-dependencies

## Call additional Makefiles in child directories
SUBDIRS = common clustalw2fasta draw-mhmm fasta-center fimo mast mcast meme mhmm mhmms motiph psp-gen \
          scaffold tomtom gomo qvalue web

DATAFILES = \
  adh.s \
  crp0.s \
  farntrans5.s \
  INO_up800.s \
  lex0.s \
  lipocalin.s \
  mini-drosoph.s \
  noq.s

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
