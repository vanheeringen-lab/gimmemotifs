## $Id: Makefile.am 4452 2010-04-14 03:45:52Z james_johnson $
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

SEDSPEC = \
	-e 's%@MEMEDIR@%$(prefix)%' \
	-e 's%@MEMEDB@%$(MEME_DB)%' \
	-e 's%@MEMELOGS@%$(MEME_LOGS)%' \
	-e 's%@SITE_URL@%$(URL)%' \
	-e 's%@METAMEME_URL@%$(URL)%' \
	-e 's%@AMIGO_URL@%$(GOSERVER_URL)%'
	
SEDBIN = \
	-e 's%@MEMEDIR@%$(prefix)%' \
	-e 's%@MEMEDB@%$(MEME_DB)%' \
	-e 's%@MEMELOGS@%$(MEME_LOGS)%' \
	-e 's%@BINARYDIR@%$(bindir)%' \
	-e 's%@APPCONFIGDIR@%$(sysconfdir)%'

SED_FIMO_OPAL = -e 's%@BINARYLOCATION@%$(FIMO_BINARY)%'
SED_GLAM2_OPAL = -e 's%@BINARYLOCATION@%$(GLAM2_BINARY)%'
SED_GLAM2SCAN_OPAL = -e 's%@BINARYLOCATION@%$(GLAM2SCAN_BINARY)%'
SED_GOMO_OPAL = -e 's%@BINARYLOCATION@%$(GOMO_BINARY)%'
SED_MAST_OPAL = -e 's%@BINARYLOCATION@%$(MAST_BINARY)%'
SED_MCAST_OPAL = -e 's%@BINARYLOCATION@%$(MCAST_BINARY)%'
SED_MEME_OPAL = -e 's%@BINARYLOCATION@%$(MEME_BINARY)%'


SHELLFILES_IN = \
  fimo.sh.in \
  glam2.sh.in \
  glam2scan.sh.in \
  gomo.sh.in \
  mast.sh.in \
  mcast.sh.in \
  meme-mast.sh.in \
	timer.sh

ETCFILES = \
    beadstring.xsl \
    bkg.jpg \
		block-diagram.xsl \
    cisml.css \
    fasta_db.csv \
    fimo.xsl \
    gomo_db.csv \
    hmm-text.xsl \
		mast-to-html.xsl \
    meme.doc \
    motif_db.csv \
    prior1.plib \
    prior30.plib \
    template.eps \
    template.svg

meme.doc: ../doc/meme.html
	elinks -dump ../doc/meme.html > meme.doc

mast.doc: ../doc/mast.html
	elinks -dump ../doc/mast.html > mast.doc

etcdir=$(sysconfdir)
install-data-hook:
if WEBSERVICE
	$(ANT) -f $(OPAL_XML) deploy -DserviceName=GOMO_${VERSION} -DappConfig=${sysconfdir}/gomo_config.xml
	$(ANT) -f $(OPAL_XML) deploy -DserviceName=FIMO_${VERSION} -DappConfig=${sysconfdir}/fimo_config.xml
	$(ANT) -f $(OPAL_XML) deploy -DserviceName=MEME_${VERSION} -DappConfig=${sysconfdir}/meme_config.xml
	$(ANT) -f $(OPAL_XML) deploy -DserviceName=MAST_${VERSION} -DappConfig=${sysconfdir}/mast_config.xml
	$(ANT) -f $(OPAL_XML) deploy -DserviceName=MCAST_${VERSION} -DappConfig=${sysconfdir}/mcast_config.xml
	$(ANT) -f $(OPAL_XML) deploy -DserviceName=GLAM2_${VERSION} -DappConfig=${sysconfdir}/glam2_config.xml
	$(ANT) -f $(OPAL_XML) deploy -DserviceName=GLAM2SCAN_${VERSION} -DappConfig=${sysconfdir}/glam2scan_config.xml
endif

uninstall-hook:
if WEBSERVICE
	$(ANT) -f $(OPAL_XML) undeploy -DserviceName=GOMO_${VERSION} -DappConfig=${sysconfdir}/gomo_config.xml
	$(ANT) -f $(OPAL_XML) undeploy -DserviceName=FIMO_${VERSION} -DappConfig=${sysconfdir}/fimo_config.xml
	$(ANT) -f $(OPAL_XML) undeploy -DserviceName=MEME_${VERSION} -DappConfig=${sysconfdir}/meme_config.xml
	$(ANT) -f $(OPAL_XML) undeploy -DserviceName=MAST_${VERSION} -DappConfig=${sysconfdir}/mast_config.xml
	$(ANT) -f $(OPAL_XML) undeploy -DserviceName=MCAST_${VERSION} -DappConfig=${sysconfdir}/mcast_config.xml
	$(ANT) -f $(OPAL_XML) undeploy -DserviceName=GLAM2_${VERSION} -DappConfig=${sysconfdir}/glam2_config.xml
	$(ANT) -f $(OPAL_XML) undeploy -DserviceName=GLAM2SCAN_${VERSION} -DappConfig=${sysconfdir}/glam2scan_config.xml
endif

etc_SCRIPTS = fimo.sh glam2.sh glam2scan.sh gomo.sh mast.sh mcast.sh meme-mast.sh timer.sh

etc_DATA = \
  $(ETCFILES) \
  cisml-to-gff.xsl \
  cisml-to-html.xsl \
  cisml-to-text.xsl \
  fimo-to-gff.xsl \
  fimo-to-html.xsl \
  fimo-to-text.xsl \
  fimo_config.xml \
  glam2_config.xml \
  glam2scan_config.xml \
  gomo_config.xml \
  gomo.xsl \
	help.gif \
	logo.js.xsl \
  mast_config.xml \
  mcast_config.xml \
  meme_config.xml \
	meme.css.xsl \
  meme.xsl

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
	
gomo.xsl: gomo.xsl.in
	$(SED) $(SEDSPEC) $? > $@
	
logo.js.xsl: logo.js
	sed -e 's/&/\&amp;/g' -e 's/</\&lt;/g' -e 's/>/\&gt;/g' $? > $?.tmp \
    && sed 's/tname/$?/g' xsl-wrapper.in | sed '3r $?.tmp' > $@ \
    && rm -f $?.tmp

mast.sh: mast.sh.in
	$(SED) $(SEDBIN) $? > $@

mast_config.xml: mast_config.xml.in
	$(SED) $(SED_MAST_OPAL) $? > $@

mcast.sh: mcast.sh.in
	$(SED) $(SEDBIN) $? > $@

mcast_config.xml: mcast_config.xml.in
	$(SED) $(SED_MCAST_OPAL) $? > $@

meme.xsl: meme.xsl.in
	$(SED) $(SEDSPEC) $? > $@

meme.css.xsl: meme.css
	sed -e 's/&/\&amp;/g' -e 's/</\&lt;/g' -e 's/>/\&gt;/g' $? > $?.tmp \
    && sed 's/tname/$?/g' xsl-wrapper.in | sed '3r $?.tmp' > $@ \
    && rm -f $?.tmp

meme_config.xml: meme_config.xml.in
	$(SED) $(SED_MEME_OPAL) $? > $@

meme-mast.sh: meme-mast.sh.in
	$(SED) $(SEDBIN) $? > $@

DATABASE_FILES = \
  db_general.csv \
  db_other_genomes.csv \
  fasta_db.csv \
  gomo_db.csv \
  motif_db.csv

# if fasta_db.csv does not exist, create an empty one.
fasta_db.csv:
	touch fasta_db.csv

# if gomo_db.csv does not exist, create an empty one.
gomo_db.csv:
	touch gomo_db.csv

# if motif_db.csv does not exist, create an empty one.
motif_db.csv:
	touch motif_db.csv

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

CLEANFILES = fimo.sh fimo_config.xml glam2.sh glam2_config.xml \
             glam2scan.sh glam2scan_config.xml gomo.xsl gomo.sh gomo_config.xml \
             mast.sh mast_config.xml mcast.sh mcast_config.xml \
             meme.xsl meme_config.xml meme-mast.sh 

EXTRA_DIST = $(SHELLFILES_IN) $(ETCFILES) $(DATABASE_FILES) \
             fimo_config.xml.in glam2_config.xml.in glam2scan_config.xml.in \
             gomo_config.xml.in gomo.xsl.in help.gif logo.js \
             mast_config.xml.in mcast_config.xml.in meme_config.xml.in meme.css meme.xsl.in \
             fimo-to-html.xsl fimo-to-text.xsl fimo-to-gff.xsl \
             cisml-to-html.xsl cisml-to-text.xsl cisml-to-gff.xsl xsl-wrapper.in


