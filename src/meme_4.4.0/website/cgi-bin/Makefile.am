## $Id: Makefile.am 4278 2009-12-23 09:58:37Z james_johnson $
##
## Process this file with automake to produce Makefile.in
##

AUTOMAKE_OPTIONS = no-dependencies

memecgidir = $(WEBDIR)/cgi-bin

SEDSPEC = -e 's%\@sendmail\@%$(SENDMAIL)%' \
          -e 's%\@mail\@%$(MAIL)%' \
	  -e 's%\@WHICHPERL\@%$(PERL)%' \
	  -e 's%\@PERLLIBDIR\@%$(libdir)/perl%' \
	  -e 's%\@contact\@%$(CONTACT)%g' \
	  -e 's%\@MEME_DIR\@%$(MEME_DIRECTORY)%' \
	  -e 's%\@MEMEDB\@%$(MEME_DB)%' \
	  -e 's%\@S_VERSION\@%$(VERSION)%' \
	  -e 's%\@SITE_URL\@%$(URL)%' \
	  -e 's%\@OPAL\@%$(OPAL_URL)%' \
	  -e 's%\@BLOCKS_URL\@%$(BLOCKS_URL)%'
	  -e 's%\@AMIGO_URL\@%$(GOSERVER_URL)%'

SUFFIXES= .pl .cgi
%.cgi: %.pl Makefile
	$(SED) $(SEDSPEC) $< > $@
%.pm: %.pm.in Makefile
	$(SED) $(SEDSPEC) $< > $@

MEME_FILES = \
  fimo.pl \
  gomo.pl \
  mast.pl \
  mcast.pl \
  meme.pl \
  get_db_list.pl \
  glam2.pl \
  glam2scan.pl \
	meme_request.pl \
  metameme.pl \
  mhmm_process_request.pl \
  nph-endpointtest.pl \
  nph-mcast-submit.pl \
  nph-placeholder.pl \
  process_request.pl \
  querystatus.pl \
  OpalServices.pm.in \
  OpalTypes.pm.in \
  submit-verify.pl \
  tomtom.pl \
  Utils.pm.in

MEME_CGI = \
  fimo.cgi \
  gomo.cgi \
  mast.cgi \
  mcast.cgi \
  meme.cgi \
	meme_request.cgi \
  get_db_list.cgi \
  glam2.cgi\
  glam2scan.cgi \
  metameme.cgi \
  mhmm_process_request.cgi  \
  nph-endpointtest.cgi \
  nph-mcast-submit.cgi \
  nph-placeholder.cgi  \
  OpalServices.pm \
  OpalTypes.pm \
  process_request.cgi \
  querystatus.cgi \
  submit-verify.cgi \
  tomtom.cgi \
  Utils.pm

if WEBSITE
memecgi_SCRIPTS = $(MEME_CGI)
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

CLEANFILES = $(MEME_CGI)

# These files should be included in the distribution.
EXTRA_DIST = $(MEME_FILES)
