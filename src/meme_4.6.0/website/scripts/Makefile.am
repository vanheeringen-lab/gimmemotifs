## $Id: Makefile.am 4911 2010-09-03 07:10:33Z james_johnson $
##
## Process this file with automake to produce Makefile.in
##

AUTOMAKE_OPTIONS = no-dependencies

moduledir=$(libdir)/perl

SEDSPEC = \
	-e 's%@SITE_CONTACT@%$(CONTACT)%'  \
	-e 's%@METAMEME_SITE_CONTACT@%$(METAMEME_CONTACT)%'  \
	-e 's%@BIN_DIR@%$(bindir)%'  \
	-e 's%@MODULE_DIR@%$(moduledir)%'  \
	-e 's%@MEMEDIR@%$(MEME_DIRECTORY)%' \
	-e 's%@METAMEMEDIR@%$(MEME_DIRECTORY)%' \
	-e 's%@MEMELOGS@%$(MEME_LOGS)%' \
	-e 's%@MEMEDB@%$(MEME_DB)%' \
	-e 's%@WEB_DIR@%$(WEBDIR)%' \
	-e 's%@METAMEME_WEB_DIR@%$(METAMEME_WEBDIR)%' \
	-e 's%@MAX_TIME@%$(MAXTIME)%' \
	-e 's%@SITE_URL@%$(URL)%' \
	-e 's%@METAMEME_SITE_URL@%$(METAMEME_URL)%'

SUFFIXES= .txt .pm
%.pm: %.txt Makefile
	$(SED) $(SEDSPEC) $< > $@

changetoweb: changetoweb.txt Makefile
	sed -e 's%#!/usr/bin/perl%#!$(PERL)%' $< > $@
metameme: metameme.txt Makefile
	sed -e 's%use lib qw(.)%use lib qw($(moduledir))%' \
      -e 's%#!/usr/bin/perl%#!$(PERL)%' $< > $@
oldmeme2meme: oldmeme2meme.txt Makefile
	sed -e 's%#!/usr/bin/perl%#!$(PERL)%' $< > $@


MEME_FILES = \
  CGIutil.pm \
  changetoweb.txt \
  Globals.txt \
  Log.txt \
  metameme.txt \
  MetaGlobals.txt \
  mhmm_globals.txt \
  oldmeme2meme.txt \
  Process.pm \
  Queue.pm \
  Submit.txt \
  UID.txt

MEME_CLEANFILES = \
  changetoweb \
  Globals.pm \
  Log.pm \
  MetaGlobals.pm \
  metameme \
  mhmm_globals.pm \
  oldmeme2meme \
  Submit.pm \
  UID.pm

if WEBSITE
MODULES_MEME = \
  Globals.pm \
  CGIutil.pm \
  Log.pm \
  MetaGlobals.pm \
  mhmm_globals.pm \
  Process.pm \
  Queue.pm \
  Submit.pm  \
  UID.pm
else 
MODULES_MEME = Globals.pm
endif

SCRIPTS_METAMEME = changetoweb metameme oldmeme2meme
module_SCRIPTS = $(MODULES_MEME) $(MODULES_METAMEME)
bin_SCRIPTS = $(SCRIPTS_METAMEME)

# The Meta-MEME web application requires the dot and convert utilities
# to be in $(bindir) due to restrictions on the PATH set by perl
# taint mode.
install-exec-hook:
if WEBSITE
	if test ! -e $(bindir)/convert; then ln -sf $(CONVERT) $(bindir)/convert; fi
	if test ! -e $(bindir)/dot; then ln -sf $(DOT) $(bindir)/dot; fi
endif

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

# These files should be deleted if the target is 'make clean'
CLEANFILES = $(MEME_CLEANFILES)

# These files should included with the distribution
EXTRA_DIST = $(MEME_FILES)
