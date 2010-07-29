## $Id: Makefile.am 4493 2010-04-22 04:38:05Z cegrant $
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
	-e 's%@SITE_URL@%$(URL)%' \
	-e 's%@DOCDIR@%$(MEME_DOC_DIR)/%'

perlmoddir=$(libdir)/perl
perlmod_SCRIPTS = convert2html.pl

MEME_TESTS = runcheck metatest

if PARALLEL
PARALLEL_TESTS = runcheck_p
if PARALLEL_IS_LAM
runcheck_p: runcheck_p.in Makefile
	sed -e "s%@mpirun@%$(top_srcdir)/scripts/lrun -n %" $< > runcheck_p; \
  chmod +x runcheck_p;
else
runcheck_p: runcheck_p.in Makefile
	sed -e "s%@mpirun@%$(MPI_CMD) %" $< > runcheck_p; \
  chmod +x runcheck_p
endif
endif

TXTFILES = \
	cat_max.txt \
	compare_dates.txt \
	fasta-fetch.txt \
	fasta-get-markov.txt \
	fasta-make-index.txt \
	fasta-shuffle-letters.txt \
	fasta-unique-names.txt \
	get_db_csv.txt \
	lrun.txt \
	mailer.txt \
	meme.txt \
	runtests \
	transfac2meme.txt \
	update-plot-usage.txt \
	update_db.txt

UTILS = \
    cat_max \
    compare_dates \
    download \
    fasta-fetch \
    fasta-get-markov \
    fasta-make-index \
    fasta-shuffle-letters \
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
    mhmm2html \
    motif2meme.pm \
    plotgen \
    priority2meme \
    prior_utils.pl \
    psp-gen \
    read_fasta_file.pl \
    sd \
    scpd2meme \
    tamo2meme \
    transfac2meme \
    uniprobe2meme \
    update_db \
    update_meme_tests

if WEBSITE
SERVER_UTILS = \
    lrun \
    mailer \
    plot-usage \
    texit \
    update-plot-usage \
    usage-reports 
else
SERVER_UTILS =
endif

# Descibe dependencies and build instructions for scripts.
# (substitute in local path to Perl).
uniprobe2meme: uniprobe2meme.pl Makefile
	$(SED) $(SEDSPEC) $< > $@
cat_max: cat_max.txt  Makefile
	$(SED) $(SEDSPEC) $< > $@
chen2meme: chen2meme.pl Makefile
	$(SED) $(SEDSPEC)  $< > $@
compare_dates: compare_dates.txt Makefile
	$(SED) $(SEDSPEC) $< > $@
convert2html.pl: convert2html.txt Makefile
	$(SED) $(SEDSPEC) $< > $@
fasta-fetch: fasta-fetch.txt Makefile
	$(SED) $(SEDSPEC) $< > $@
fasta-get-markov: fasta-get-markov.txt Makefile
	$(SED) $(SEDSPEC) $< > $@
fasta-make-index: fasta-make-index.txt Makefile
	$(SED) $(SEDSPEC) $< > $@
fasta-shuffle-letters: fasta-shuffle-letters.txt Makefile
	$(SED) $(SEDSPEC) $< > $@
fasta-unique-names: fasta-unique-names.txt Makefile
	$(SED) $(SEDSPEC) $< > $@
get_db_csv: get_db_csv.txt Makefile
	$(SED) $(SEDSPEC) $< > $@
glam2html: glam2html.py Makefile
	$(SED) $(SEDSPEC) $< > $@
glam2scan2html: glam2scan2html.py Makefile
	$(SED) $(SEDSPEC) $< > $@
glam2psfm: glam2psfm.py Makefile
	$(SED) $(SEDSPEC)  $< > $@
hart2meme-bkg: hart2meme-bkg.pl Makefile
	$(SED) $(SEDSPEC)  $< > $@
hartemink2psp: hartemink2psp.pl Makefile
	$(SED) $(SEDSPEC) $< > $@
iupac2meme: iupac2meme.pl Makefile
	$(SED) $(SEDSPEC) $< > $@
jaspar2meme: jaspar2meme.pl Makefile
	$(SED) $(SEDSPEC) $< > $@
lrun: lrun.txt Makefile
	$(SED) $(SEDSPEC) $< > $@; \
  chmod +x $@
psp-gen: psp-gen.pl prior_utils.pl \
 read_fasta_file.pl Makefile
	$(SED) $(SEDSPEC) $< > $@; \
  chmod +x $@
mailer: mailer.txt Makefile
	$(SED) $(SEDSPEC) $< > $@
meme: meme.txt Makefile
	$(SED) $(SEDSPEC) $< > $@
metatest: metatest.pl Makefile
	$(SED) $(SEDSPEC) $< > $@
	chmod ugo+x metatest
mhmm2html: mhmm2html.pl Makefile
	$(SED) $(SEDSPEC) $< > $@
motif2meme.pm: motif2meme.pm.in Makefile
	$(SED) $(SEDSPEC) $< > $@
priority2meme: priority2meme.pl
	$(SED) $(SEDSPEC) $? > $@
scpd2meme: scpd2meme.pl Makefile
	$(SED) $(SEDSPEC) $< > $@
tamo2meme: tamo2meme.pl Makefile
	$(SED) $(SEDSPEC) $< > $@
transfac2meme: transfac2meme.txt Makefile
	$(SED) $(SEDSPEC) $< > $@
update_db: update_db.txt  Makefile
	$(SED) $(SEDSPEC) $< > $@
update-plot-usage: update-plot-usage.txt Makefile
	$(SED) $(SEDSPEC) $<  > $@

# Set automake targets
TESTS_ENVIRONMENT = testdir=$(top_srcdir)/tests
TESTS = $(MEME_TESTS) $(PARALLEL_TESTS)
bin_SCRIPTS = $(UTILS) $(SERVER_UTILS)

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
  uniprobe \
  convert2html.pl \
  fasta-get-markov \
  glam2html \
  glam2psfm \
  glam2scan2html \
  hart2meme \
  hartemink2psp \
  iupac2meme \
  jaspar2meme \
  psp-gen \
  iupac2meme.pm \
  priority2meme \
  scpd2meme \
  tamo2meme \
  runcheck_p

EXTRA_DIST = \
  uniprobe2meme.pl \
  chen2meme.pl \
  convert2html.txt \
  download \
  glam2html.py \
  glam2psfm.py \
  glam2scan2html.py \
  hartemink2psp.pl \
  hart2meme-bkg.pl \
  iupac2meme.pl \
  jaspar2meme.pl \
  psp-gen.pl \
  metatest.pl \
  mhmm2html.pl \
  motif2meme.pm.in \
	obo2dag.jar \
  plotgen \
  plot-usage \
  priority2meme.pl \
  prior_utils.pl \
  read_fasta_file.pl \
  scpd2meme.pl \
  sd \
  tamo2meme.pl \
  $(TXTFILES) \
  runcheck \
  runcheck_p.in \
  texit \
  update_meme_tests \
  usage-reports

dist-hook:
	rm -f $(distdir)/metatest
	rm -f $(distdir)/runcheck_p
