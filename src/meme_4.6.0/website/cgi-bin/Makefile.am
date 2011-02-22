## $Id: Makefile.am 5298 2011-01-17 05:48:02Z cegrant $
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
  meme-chip.pl \
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
  spamo.pl \
  tomtom.pl \
  tomtom_request.pl

WEBSITE_TEMPLATES = \
  db_list_view.tmpl \
  motif_input.tmpl \
  spamo.tmpl \
  spamo_verify.tmpl \
  tomtom.tmpl \
  tomtom_verify.tmpl

# Templates on both sides of a split install
DUALUSE_TEMPLATES = \
  job_status.tmpl

WEBSERVICE_TEMPLATES = 

WEBSITE_CGI = \
  fimo.cgi \
  gomo.cgi \
  mast.cgi \
  mcast.cgi \
  meme.cgi \
  meme-chip.cgi \
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
  spamo.cgi \
  submit-verify.cgi \
  tomtom.cgi \
  tomtom_request.cgi

MEME_TEMPLATES =
MEME_CGI =

# WEBSITE OR WEBSERVICE
if WEBSITE
MEME_TEMPLATES += $(DUALUSE_TEMPLATES)
else !WEBSITE
if WEBSERVICE
MEME_TEMPLATES += $(DUALUSE_TEMPLATES)
endif WEBSERVICE
endif !WEBSITE

if WEBSITE
MEME_TEMPLATES += $(WEBSITE_TEMPLATES)
MEME_CGI += $(WEBSITE_CGI)
endif

if WEBSERVICE
MEME_TEMPLATES += $(WEBSERVICE_TEMPLATES)
endif

memecgi_SCRIPTS = $(MEME_CGI)
memecgi_DATA = $(MEME_TEMPLATES)

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
EXTRA_DIST = $(MEME_FILES) $(WEBSITE_TEMPLATES) $(DUALUSE_TEMPLATES) $(WEBSERVICE_TEMPLATES)
