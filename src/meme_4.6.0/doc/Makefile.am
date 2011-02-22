#####################################################################
## FILE: Makefile.am
# AUTHOR: William Stafford Noble, Charles E. Grant
# CREATE DATE: 2003
# PROJECT: Meta-MEME
# COPYRIGHT: 2003, UW
#####################################################################

## Automake uses this file to generate a Makefile

#####################################################################
#  This is the Makefile for the documentation directory.
#####################################################################

SUBDIRS = examples images .

SEDSPEC = -e 's%\@VERSION\@%$(VERSION)%g' \
          -e 's%\@CONTACT\@%$(CONTACT)%g' \
          -e 's%\@METAMEME_CONTACT\@%$(METAMEME_CONTACT)%g' \
          -e 's%\@URL\@%$(URL)%g'

# Rule for creating HTML output from XML input.
# Uses XSLT to expand common option descriptions.
# Can also be used to add other common elements.
%.html: %.xml options.xsl
	../src/xsltproc_lite $(srcdir)/options.xsl $< $@

modified-template-footer.js: modified-template-footer.in Makefile
	$(SED) $(SEDSPEC) $< > $@

# XML source files
# Included in distribution, but not installed.
COMMON_XML = \
  alphtype.xml \
  ama.xml \
  ame.xml \
  beadstring.xml \
  ceqlogo.xml \
  clustalw2phylip.xml \
  compute-prior-dist.xml \
  draw-mhmm.xml \
  dreme.xml \
  fasta-get-markov.xml \
  fasta-io.xml \
  fimo.xml \
  fitevd.xml \
  gendb.xml \
  getsize.xml \
  gomo.xml \
  log-hmm.xml \
  mcast.xml \
  meme-chip-tutorial.xml \
  meme-io.xml \
  mhmm.xml \
  mhmme.xml \
  mhmms.xml \
  miao.xml \
  motiph.xml \
  options.xsl \
  pmp_bf.xml \
  psp-gen.xml \
  qvalue.xml \
  readseq.xml \
  reconcile-tree-alignment.xml \
  reduce-alignment.xml \
  remove-alignment-gaps.xml \
  shadow.xml

# Included in distribution, but not installed.
COMMON_IN_FILES = \
  modified-template-footer.in

# Files generated from XML source
# Not included in distribution, but installed.
COMMON_BUILT_DOC = \
  alphtype.html \
  ama.html \
  ame.html \
  beadstring.html \
  ceqlogo.html \
  clustalw2phylip.html \
  compute-prior-dist.html \
  draw-mhmm.html \
  dreme.html \
  fasta-get-markov.html \
  fasta-io.html \
  fimo.html \
  fitevd.html \
  gendb.html \
  getsize.html \
  gomo.html \
  log-hmm.html \
  mcast.html \
  meme-chip-tutorial.html \
  meme-io.html \
  mhmm.html \
  mhmme.html \
  mhmms.html \
  miao.html \
  modified-template-footer.js \
  motiph.html \
  pmp_bf.html \
  psp-gen.html \
  qvalue.html \
  readseq.html \
  reconcile-tree-alignment.html \
  reduce-alignment.html \
  remove-alignment-gaps.html \
  shadow.html

# Files Included in distribution and installed.
COMMON_DOC = \
  alphabets.html \
  authors.html \
  bfile-format.html \
  chen2meme.html \
  cite.html \
  clustalw2fasta.html \
  clustalw-format.html \
  copyright.html \
  fasta-format.html \
  glam2.css \
  glam2format_man.html \
  glam2mask_man.html \
  glam2scan_man.html \
  glam2_alphabet.html \
  glam2_dirichlet.html \
  glam2_install.html \
  glam2_man.html \
  glam2-template-header.js \
  GLAM2_method.pdf \
  glam2_ref.html \
  glam2_tut.html \
  godag-format.html \
  gomo-format.html \
  gomo_highlight.html \
  install.html \
  iupac.html \
  iupac2meme.html \
  jaspar2meme.html \
  lipo.fasta \
  lipo.linear.all.mhmms \
  lipo.linear.mhmm \
  lipo.meme.html \
  lipo.meme.stripped.txt \
  mast.html \
	mast2txt.html \
  meme.dtd \
  meme.html \
  meme-format.html \
  meme-install.html \
  meme-suite-logo.js \
  meme-suite.css \
  metameme.css \
  meta-meme.gif \
  mhmm-output.html \
  mhmmscan.html \
  overview.html \
  priority2meme.html \
  purge.html \
  release-notes.html \
  scpd2meme.html \
  scripting.html \
  sdsc.gif \
  select-star-species.html \
  standards.html \
  tamo2meme.html \
  tomtom.html \
  transfac2meme.html \
  transfac-format.html \
  uniprobe2meme.html

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


# Set installation location for doc files.
if WEBSITE
MEME_DOC_DIR = $(WEBDIR)/doc
else
MEME_DOC_DIR = $(prefix)/doc
endif

memedocdir = $(MEME_DOC_DIR)
memedoc_DATA =  $(COMMON_DOC) $(COMMON_BUILT_DOC) $(MEME_DOC) $(MEME_BUILT_DOC)

CLEANFILES = $(COMMON_BUILT_DOC) 
EXTRA_DIST = $(COMMON_XML) $(COMMON_IN_FILES) $(COMMON_DOC)
