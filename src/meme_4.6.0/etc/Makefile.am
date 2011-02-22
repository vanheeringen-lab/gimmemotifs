## $Id: Makefile.am 5305 2011-01-20 21:52:46Z cegrant $
##
## Process this file with automake to produce Makefile.in
##

AUTOMAKE_OPTIONS = no-dependencies

MEME_BINARY_DIR = $(prefix)/bin
FIMO_BINARY = $(sysconfdir)/fimo.sh
GLAM2_BINARY = $(sysconfdir)/glam2.sh
GLAM2SCAN_BINARY = $(sysconfdir)/glam2scan.sh
GOMO_BINARY = $(sysconfdir)/gomo.sh
MAST_BINARY = $(sysconfdir)/mast.sh
MCAST_BINARY = $(sysconfdir)/mcast.sh
MEME_BINARY = $(sysconfdir)/meme-mast.sh
MEMECHIP_BINARY = $(sysconfdir)/meme-chip
TOMTOM_SERVICE = $(bindir)/tomtom_webservice
SPAMO_SERVICE = $(bindir)/spamo_webservice

SEDSPEC = \
	-e 's%@MEMEDIR@%$(prefix)%' \
	-e 's%@MEMEDB@%$(MEME_DB)%' \
	-e 's%@MEMELOGS@%$(MEME_LOGS)%' \
	-e 's%@SITE_URL@%$(URL)%' \
	-e 's%@METAMEME_URL@%$(URL)%' \
	-e 's%@AMIGO_URL@%$(GOSERVER_URL)%'
	
SEDBIN = \
	-e 's%@WHICHPERL@%$(PERL)%' \
	-e 's%@PERLMODDIR@%$(perlmoddir)%' \
	-e 's%@PERLLIBDIR@%$(libdir)/perl%' \
	-e 's%@MEMEDIR@%$(prefix)%' \
	-e 's%@MEMEDB@%$(MEME_DB)%' \
	-e 's%@MEMELOGS@%$(MEME_LOGS)%' \
	-e 's%@BINARYDIR@%$(bindir)%' \
	-e 's%@APPCONFIGDIR@%$(sysconfdir)%' \
	-e 's%@SITE_URL@%$(URL)%'

SED_FIMO_OPAL = -e 's%@BINARYLOCATION@%$(FIMO_BINARY)%'
SED_GLAM2_OPAL = -e 's%@BINARYLOCATION@%$(GLAM2_BINARY)%'
SED_GLAM2SCAN_OPAL = -e 's%@BINARYLOCATION@%$(GLAM2SCAN_BINARY)%'
SED_GOMO_OPAL = -e 's%@BINARYLOCATION@%$(GOMO_BINARY)%'
SED_MEME_OPAL = -e 's%@BINARYLOCATION@%$(MEME_BINARY)%'
SED_MEMECHIP_OPAL = -e 's%@BINARYLOCATION@%$(MEMECHIP_BINARY)%'
SED_MAST_OPAL = -e 's%@BINARYLOCATION@%$(MAST_BINARY)%'
SED_MCAST_OPAL = -e 's%@BINARYLOCATION@%$(MCAST_BINARY)%'
SED_TOMTOM_OPAL = -e 's%@TOMTOM_SERVICE@%$(TOMTOM_SERVICE)%'
SED_SPAMO_OPAL = -e 's%@SPAMO_SERVICE@%$(SPAMO_SERVICE)%'


SHELLFILES_IN = \
  fimo.sh.in \
  glam2.sh.in \
  glam2scan.sh.in \
  gomo.sh.in \
  mast.sh.in \
  mcast.sh.in \
  meme-mast.sh.in \
  meme-chip.txt \
  timer.sh

ETCFILES = \
  beadstring.xsl \
  bkg.jpg \
  block-diagram.xsl \
  cisml.css \
  cisml-to-gff3.xsl \
  cisml-to-html.xsl \
  cisml-to-text.xsl \
  fimo-to-wiggle.xsl \
  fimo-to-gff3.xsl \
  fimo-to-html.xsl \
  fimo-to-text.xsl \
  gomo-to-html.xsl \
  help.gif \
  hmm-text.xsl \
  mast-to-html.xsl \
  meme-to-html.xsl \
  meme.doc \
  prior1.plib \
  prior30.plib \
  spamo_template.eps \
  spamo-to-html.xsl \
  template.eps \
  template.svg \
  tomtom-to-html.xsl

meme.doc: ../doc/meme.html
	elinks -dump ../doc/meme.html > meme.doc

mast.doc: ../doc/mast.html
	elinks -dump ../doc/mast.html > mast.doc

etcdir=$(sysconfdir)
install-data-hook:
if WEBSERVICE
	$(ANT) -f $(OPAL_XML) deploy -DserviceName=FIMO_${VERSION} -DappConfig=${sysconfdir}/fimo_config.xml
	$(ANT) -f $(OPAL_XML) deploy -DserviceName=GOMO_${VERSION} -DappConfig=${sysconfdir}/gomo_config.xml
	$(ANT) -f $(OPAL_XML) deploy -DserviceName=MAST_${VERSION} -DappConfig=${sysconfdir}/mast_config.xml
	$(ANT) -f $(OPAL_XML) deploy -DserviceName=MCAST_${VERSION} -DappConfig=${sysconfdir}/mcast_config.xml
	$(ANT) -f $(OPAL_XML) deploy -DserviceName=MEME_${VERSION} -DappConfig=${sysconfdir}/meme_config.xml
	$(ANT) -f $(OPAL_XML) deploy -DserviceName=MEMECHIP_${VERSION} -DappConfig=${sysconfdir}/meme-chip_config.xml
	$(ANT) -f $(OPAL_XML) deploy -DserviceName=GLAM2_${VERSION} -DappConfig=${sysconfdir}/glam2_config.xml
	$(ANT) -f $(OPAL_XML) deploy -DserviceName=GLAM2SCAN_${VERSION} -DappConfig=${sysconfdir}/glam2scan_config.xml
	$(ANT) -f $(OPAL_XML) deploy -DserviceName=TOMTOM_${VERSION} -DappConfig=${sysconfdir}/tomtom_config.xml
	$(ANT) -f $(OPAL_XML) deploy -DserviceName=SPAMO_${VERSION} -DappConfig=${sysconfdir}/spamo_config.xml
endif

uninstall-hook:
if WEBSERVICE
	$(ANT) -f $(OPAL_XML) undeploy -DserviceName=FIMO_${VERSION} -DappConfig=${sysconfdir}/fimo_config.xml
	$(ANT) -f $(OPAL_XML) undeploy -DserviceName=GOMO_${VERSION} -DappConfig=${sysconfdir}/gomo_config.xml
	$(ANT) -f $(OPAL_XML) undeploy -DserviceName=MAST_${VERSION} -DappConfig=${sysconfdir}/mast_config.xml
	$(ANT) -f $(OPAL_XML) undeploy -DserviceName=MCAST_${VERSION} -DappConfig=${sysconfdir}/mcast_config.xml
	$(ANT) -f $(OPAL_XML) undeploy -DserviceName=MEME_${VERSION} -DappConfig=${sysconfdir}/meme_config.xml
	$(ANT) -f $(OPAL_XML) undeploy -DserviceName=MEMECHIP_${VERSION} -DappConfig=${sysconfdir}/meme-chip_config.xml
	$(ANT) -f $(OPAL_XML) undeploy -DserviceName=GLAM2_${VERSION} -DappConfig=${sysconfdir}/glam2_config.xml
	$(ANT) -f $(OPAL_XML) undeploy -DserviceName=GLAM2SCAN_${VERSION} -DappConfig=${sysconfdir}/glam2scan_config.xml
	$(ANT) -f $(OPAL_XML) undeploy -DserviceName=TOMTOM_${VERSION} -DappConfig=${sysconfdir}/tomtom_config.xml
	$(ANT) -f $(OPAL_XML) undeploy -DserviceName=SPAMO_${VERSION} -DappConfig=${sysconfdir}/spamo_config.xml
endif

etc_SCRIPTS = fimo.sh glam2.sh glam2scan.sh gomo.sh mast.sh mcast.sh meme-mast.sh meme-chip timer.sh

etc_DATA = \
  $(DATABASE_FILES) \
  $(ETCFILES) \
  constants.xsl \
  delay_draw.js.xsl \
  fimo_config.xml \
  glam2_config.xml \
  glam2scan_config.xml \
  gomo_config.xml \
  motif_logo.js.xsl \
  mast_config.xml \
  mcast_config.xml \
  meme_config.xml \
  meme-chip_config.xml \
  meme.css.xsl \
  spamo-to-html.css.xsl \
  spamo-to-html.js.xsl \
  spamo_graph.js.xsl \
  spamo_config.xml \
  tomtom_config.xml

if WEBSITE
memehtmldir = $(WEBDIR)
# _DATA doesn't automaticaly include the file in the distribution so it needs to be in EXTRA_DIST too
memehtml_DATA = \
  motif_logo.js \
  meme.css
endif

constants.xsl: constants.xsl.in
	$(SED) $(SEDSPEC) $? > $@

%.xsl: %
	$(SED) 's/tname/$?/' xsl-wrapper.in | $(SED) '5r $?' > $@

fimo.sh: fimo.sh.in
	$(SED) $(SEDBIN) $? > $@

fimo_config.xml: fimo_config.xml.in
	$(SED) $(SED_FIMO_OPAL) $? > $@

glam2.sh: glam2.sh.in
	$(SED) $(SEDBIN) $? > $@

glam2_config.xml: glam2_config.xml.in
	$(SED) $(SED_GLAM2_OPAL) $? > $@

glam2scan_config.xml: glam2scan_config.xml.in
	$(SED) $(SED_GLAM2SCAN_OPAL) $? > $@

glam2scan.sh: glam2scan.sh.in
	$(SED) $(SEDBIN) $? > $@

gomo.sh: gomo.sh.in
	$(SED) $(SEDBIN) $? > $@

gomo_config.xml: gomo_config.xml.in
	$(SED) $(SED_GOMO_OPAL) $? > $@

mast.sh: mast.sh.in
	$(SED) $(SEDBIN) $? > $@

mast_config.xml: mast_config.xml.in
	$(SED) $(SED_MAST_OPAL) $? > $@

mcast.sh: mcast.sh.in
	$(SED) $(SEDBIN) $? > $@

mcast_config.xml: mcast_config.xml.in
	$(SED) $(SED_MCAST_OPAL) $? > $@

meme_config.xml: meme_config.xml.in
	$(SED) $(SED_MEME_OPAL) $? > $@

meme-chip_config.xml: meme-chip_config.xml.in
	$(SED) $(SED_MEMECHIP_OPAL) $? > $@

meme-mast.sh: meme-mast.sh.in
	$(SED) $(SEDBIN) $? > $@

meme-chip: meme-chip.txt
	$(SED) $(SEDBIN) $? > $@

spamo_config.xml: spamo_config.xml.in
	$(SED) $(SED_SPAMO_OPAL) $? > $@

tomtom_config.xml: tomtom_config.xml.in
	$(SED) $(SED_TOMTOM_OPAL) $? > $@

DATABASE_INDEXES = \
  fasta_db.index \
  gomo_db.index \
  motif_db.index

DATABASE_FILES = \
  db_general.csv \
  db_other_genomes.csv \
  fasta_db.csv \
  gomo_db.csv \
  motif_db.csv \
  $(DATABASE_INDEXES)

# if database does not exist, create an empty one.
%_db.csv:
	touch $@

# index the database
%_db.index: %_db.csv
	$(PERL) -I ../scripts/ -MCatList -e 'CatList::index_catlist_csv("$?", "$@")'


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

CLEANFILES = constants.xsl delay_draw.js.xsl fimo.sh fimo_config.xml glam2.sh \
             glam2_config.xml glam2scan.sh glam2scan_config.xml gomo.sh \
             gomo_config.xml gomo_db.index mast.sh mast_config.xml mcast.sh \
             mcast_config.xml meme.css.xsl meme_config.xml meme-chip_config.xml \
             meme-chip  meme-mast.sh motif_logo.js.xsl spamo-to-html.css.xsl \
             spamo-to-html.js.xsl spamo_graph.js.xsl $(DATABASE_INDEXES)

EXTRA_DIST = $(SHELLFILES_IN) $(ETCFILES) \
             constants.xsl.in db_general.csv db_other_genomes.csv delay_draw.js \
             fimo_config.xml.in glam2_config.xml.in glam2scan_config.xml.in \
             gomo_config.xml.in motif_logo.js mast_config.xml.in \
             mcast_config.xml.in meme_config.xml.in  meme-chip_config.xml.in \
             meme.css spamo-to-html.css spamo-to-html.js spamo_graph.js \
             spamo_config.xml.in tomtom_config.xml.in \
             xsl-wrapper.in


