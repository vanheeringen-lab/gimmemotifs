## $Id: Makefile.am 3890 2009-07-18 09:28:32Z tbailey $
##
## Process this file with automake to produce Makefile.in

AUTOMAKE_OPTIONS = no-dependencies

SUBDIRS = src scripts tests website doc etc

install-data-hook: dbdir

dbdir:
	mkdir -p $(MEME_DB)

install-data-local:
	mkdir -p $(MEME_LOGS)
	chmod a+w $(MEME_LOGS)

distdir = $(PACKAGE)_$(VERSION)

distclean-local: clean-local
clean-local: clean-ac
clean-ac:
	rm -rf autom4te.cache

test:
	cd scripts; $(MAKE) check

# Do not include .svn directories in the distribution
dist-hook:
	rm -rf `find $(distdir) -name .svn`
