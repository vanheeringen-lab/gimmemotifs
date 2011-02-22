## $Id: Makefile.am 5280 2011-01-05 11:10:34Z james_johnson $
##
## Process this file with automake to produce Makefile.in
##
##

AUTOMAKE_OPTIONS = no-dependencies

SEDSPEC = -e 's%\@VERSION\@%$(VERSION)%g' \
          -e 's%\@CONTACT\@%$(CONTACT)%g' \
          -e 's%\@METAMEME_CONTACT\@%$(METAMEME_CONTACT)%g' \
          -e 's%\@URL\@%$(URL)%g'

SUFFIXES = html
HEAD = head-copyright
TAIL = tail-copyright

COMMON_INSTALL = $(COMMON_BUILT_FILES)

memehtmldir = $(WEBDIR)

MEME_FILES = \
  alternates.html \
  canvas-test.html \
  check-submission-form.js \
  cisml.css \
  cisml.pattern.html.xsl \
  documentation.html \
  downloads.html \
  fimo-intro.html \
  general-faq.html \
  glam2-intro.html \
  glam2scan-intro.html \
  gomo-intro.html \
  $(HEAD) \
  help_alphabet.html \
  help_distribution.html \
  help_email.html \
  help_filters.html \
  help_format.html \
  help_multiple_genomes.html \
  help_psp.html \
  help_sequences.html \
  help_width.html \
  intro.html \
  mast-faq.html \
  mast-input.html \
  mast-intro.html \
  mast-output.html \
  mcast-intro.html \
  meme-download.html \
  meme-faq.html \
  meme-input.html \
  memechip-intro.html \
  meme-intro.html \
  metameme-help-format.html \
  metameme-intro.html \
  metameme.css \
  motif_discovery.html \
	motif_input.css \
  motif_input.js \
  motif_search.html \
  papers.html \
  sample.database \
  sample.query \
  search.js \
  selectdb.js \
  showhide.js \
  site.css \
	spamo.js \
  $(TAIL) \
  template-css.js \
  template-header.js \
  tomtom-intro.html \
  tomtom.js \
  validate.js

MEME_IN_FILES = \
  cisml.pattern.html.xsl.in \
  meme-suite-menu.in \
  metameme-help-format.in \
  metameme-intro.in \
  mhmm-download.in \
  mhmm-faq.in \
  mhmm-links.in \
  mhmm-options.in \
  mhmm-overview.in \
  mhmm-papers.in \
  resources.in \
  search.in \
  submit-verify.in \
  template-footer.in

# These files must be built from the corresponding .in files
# and should not be included in the distribution.
MEME_BUILT_FILES = \
  cisml.pattern.html.xsl \
  COPYRIGHT.html \
  meme-suite-menu.js \
  metameme-help-format.html \
  metameme-intro.html \
  mhmm-download.html \
  mhmm-faq.html \
  mhmm-links.html \
  mhmm-options.html \
  mhmm-overview.html \
  mhmm-papers.html \
  submit-verify.html \
  resources.html \
  search.html \
  template-footer.js

cisml.pattern.html.xsl: cisml.pattern.html.xsl.in Makefile
	$(SED) $(SEDSPEC) cisml.pattern.html.xsl.in > $@

COPYRIGHT.html: $(HEAD) $(top_srcdir)/COPYING
	/bin/cat $(HEAD) $(top_srcdir)/COPYING $(TAIL) > COPYRIGHT.html

meme-suite-menu.js: meme-suite-menu.in Makefile
	$(SED)  $(SEDSPEC) meme-suite-menu.in > $@

metameme-help-format.html: metameme-help-format.in Makefile
	$(SED) $(SEDSPEC) metameme-help-format.in > $@

metameme-intro.html: metameme-intro.in Makefile
	$(SED) $(SEDSPEC) metameme-intro.in > $@

mhmm-download.html: mhmm-download.in Makefile
	$(SED) $(SEDSPEC) mhmm-download.in > $@

mhmm-faq.html: mhmm-faq.in Makefile
	$(SED) $(SEDSPEC) mhmm-faq.in > $@

mhmm-links.html: mhmm-links.in Makefile
	$(SED) $(SEDSPEC) mhmm-links.in > $@

mhmm-options.html: mhmm-options.in Makefile
	$(SED) $(SEDSPEC) mhmm-options.in > $@

mhmm-overview.html: mhmm-overview.in Makefile
	$(SED) $(SEDSPEC) mhmm-overview.in > $@

mhmm-papers.html: mhmm-papers.in Makefile
	$(SED) $(SEDSPEC) mhmm-papers.in > $@

resources.html: resources.in Makefile
	$(SED) $(SEDSPEC) resources.in > $@

search.html: search.in Makefile
	$(SED) $(SEDSPEC) search.in > $@

submit-verify.html: submit-verify.in Makefile
	$(SED) $(SEDSPEC) submit-verify.in > $@

template-footer.js: template-footer.in Makefile
	$(SED)  $(SEDSPEC) template-footer.in > $@

MEME_INSTALL = $(MEME_FILES) $(MEME_BUILT_FILES) COPYRIGHT.html


if WEBSITE

install-data-hook:
	cd $(WEBDIR); \
	if [ ! -e index.html ]; then \
		ln -s intro.html index.html; \
	fi

endif

memehtml_DATA = $(MEME_FILES) $(MEME_BUILT_FILES)

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

CLEANFILES = \
  $(MEME_BUILT_FILES) \
  index.html

EXTRA_DIST = \
  $(MEME_FILES) \
  $(MEME_IN_FILES) \
  $(HEAD) \
  $(TAIL)
