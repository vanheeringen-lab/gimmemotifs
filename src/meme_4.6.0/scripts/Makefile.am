## $Id: Makefile.am 5294 2011-01-11 03:07:15Z james_johnson $
##
## Process this file with automake to produce Makefile.in
##
##

# Set installation location for doc files.
if WEBSITE
MEME_DOC_DIR = $(WEBDIR)/doc
else
MEME_DOC_DIR = $(prefix)/doc
endif

SEDSPEC = \
	-e 's%@WEB_DIR@%$(WEBDIR)%' \
	-e 's%@WHICHPERL@%$(PERL)%' \
	-e 's%@PERLMODDIR@%$(perlmoddir)%' \
	-e 's%@PERLLIBDIR@%$(libdir)/perl%' \
	-e 's%@WHICHPYTHON@%$(PYTHON)%' \
	-e 's%@MACH_ID@%$(MACHID)%' \
	-e 's%@sendmail@%$(SENDMAIL)%' \
	-e 's%@mail@%$(MAIL)%' \
	-e 's%@CONTACT@%$(CONTACT)%' \
	-e 's%@LAM_BIN@%$(LAMBIN)%' \
	-e 's%@MPICMD@%$(MPI_CMD)%' \
	-e 's%@MEMEDIR@%$(MEME_DIRECTORY)%' \
	-e 's%@MEME_ETC_DIR@%$(top_srcdir)/etc%' \
	-e 's%@BINDIR@%$(bindir)%' \
	-e 's%@MEMELOGS@%$(MEME_LOGS)%' \
	-e 's%@MEMEDB@%$(MEME_DB)%' \
	-e 's%@HOST_NAME@%$(HOSTNAME)%' \
	-e 's%@MEMEUSER@%$(MEME_USER)%' \
	-e 's%@MPI_TYPE@%$(MPITYPE)%' \
	-e 's%@S_VERSION@%$(VERSION)%' \
	-e 's%@SITE_URL@%$(URL)%' \
	-e 's%@OPAL@%$(OPAL_URL)%' \
	-e 's%@DOCDIR@%$(MEME_DOC_DIR)/%'

BUILT_PERL_MOD = MemeWebUtils.pm convert2html.pl

perlmoddir=$(libdir)/perl
perlmod_SCRIPTS = $(BUILT_PERL_MOD) MotifUtils.pm CatList.pm

MEME_TESTS = runcheck metatest

if PARALLEL
PARALLEL_TESTS = runcheck_p
if PARALLEL_IS_LAM
runcheck_p: runcheck_p.sh.in Makefile
	sed -e "s%@mpirun@%$(top_srcdir)/scripts/lrun -n %" $< > runcheck_p; \
  chmod +x runcheck_p;
else
runcheck_p: runcheck_p.sh.in Makefile
	sed -e "s%@mpirun@%$(MPI_CMD) %" $< > runcheck_p; \
  chmod +x runcheck_p
endif
endif


INPUT_FILES = $(wildcard *.in)

BUILT_UTILS = \
  ama-qvalues \
  cat_max \
  compare_dates \
  dreme \
  fasta-center \
  fasta-dinucleotide-shuffle \
  fasta-dinucleotide-shuffle.py \
  fasta-fetch \
  fasta-get-markov \
  fasta-hamming-enrich \
  fasta-hamming-enrich.py \
  fasta-make-index \
  fasta-shuffle-letters \
  fasta-subsample \
  fasta-unique-names \
  glam2html \
  glam2psfm \
  glam2scan2html \
  get_db_csv \
  chen2meme \
  jaspar2meme \
  hart2meme-bkg \
  hartemink2psp \
  iupac2meme \
  meme \
  meme-get-motif \
  meme2meme \
  mhmm2html \
  nmica2meme \
  priority2meme \
  psp-gen \
  rna2meme \
  scpd2meme \
  taipale2meme \
  tamo2meme \
  transfac2meme \
  uniprobe2meme \
  update_db

UTILS = \
  download \
  hypergeometric.py \
  plotgen \
  prior_utils.pl \
  read_fasta_file.pl \
  sd \
  sequence.py \
  update_meme_tests \
  $(BUILT_UTILS)

# services
DUALUSE_SERVICES = \
  tomtom_webservice

WEBSERVICE_SERVICES = \
  spamo_webservice

#utilities
WEBSITE_BUILT_WEBUTILS = \
  lrun \
  mailer \
  update-plot-usage

WEBSITE_WEBUTILS = \
  plot-usage \
  texit \
  usage-reports \
  $(WEBSITE_BUILT_WEBUTILS)

MEME_SERVICES =
MEME_WEBUTILS =

# WEBSITE OR WEBSERVICE
if WEBSITE
MEME_SERVICES += $(DUALUSE_SERVICES)
else !WEBSITE
if WEBSERVICE
MEME_SERVICES += $(DUALUSE_SERVICES)
endif WEBSERVICE
endif !WEBSITE


if WEBSERVICE
MEME_SERVICES += $(WEBSERVICE_SERVICES)
endif

if WEBSITE
MEME_WEBUTILS += $(WEBSITE_WEBUTILS)
endif


# Describe dependencies and build instructions for scripts.
# (substitute in local path to Perl).

# input files with unusual build names
convert2html.pl: convert2html.pl.in
	$(SED) $(SEDSPEC) $< > $@; \
	chmod +x $@
MemeWebUtils.pm: MemeWebUtils.pm.in Makefile
	$(SED) $(SEDSPEC) $< > $@

# need these to allow use as a Python library or for command line
fasta-hamming-enrich.py: fasta-hamming-enrich
	cp fasta-hamming-enrich fasta-hamming-enrich.py
fasta-dinucleotide-shuffle.py: fasta-dinucleotide-shuffle
	cp fasta-dinucleotide-shuffle fasta-dinucleotide-shuffle.py

# csh input files
%: %.csh.in
	$(SED) $(SEDSPEC) $< > $@; \
  chmod +x $@

# perl input files
%: %.pl.in
	$(SED) $(SEDSPEC) $< > $@; \
  chmod +x $@

# python input files
%: %.py.in sequence.py hypergeometric.py
	$(SED) $(SEDSPEC) $< > $@; \
  chmod +x $@


# Set automake targets
TESTS_ENVIRONMENT = testdir=$(top_srcdir)/tests
TESTS = $(MEME_TESTS) $(PARALLEL_TESTS)
bin_SCRIPTS = $(UTILS) $(MEME_WEBUTILS) $(MEME_SERVICES)

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

svndir=$(prefix)

# Define which files to delete if target is 'make clean'.
CLEANFILES = \
  $(BUILT_PERL_MOD) \
  $(BUILT_UTILS) \
  $(WEBSITE_BUILT_WEBUTILS) \
  $(MEME_SERVICES) \
  metatest

EXTRA_DIST = \
  CatList.pm \
  comb-dynamicprog.pl \
  download \
  hypergeometric.py \
  hypergeometric-dynamicprog.pl \
  MotifUtils.pm \
  obo2dag.jar \
  plotgen \
  plot-usage \
  prior_utils.pl \
  read_fasta_file.pl \
  runcheck \
  runtests \
  sd \
  sequence.py \
  texit \
  update_meme_tests \
  usage-reports \
  $(INPUT_FILES)

dist-hook:
	rm -f $(distdir)/metatest
	rm -f $(distdir)/runcheck_p
