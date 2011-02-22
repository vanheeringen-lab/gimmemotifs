## $Id: Makefile.am 2599 2008-04-30 08:59:29Z tbailey $
##
## Process this file with automake to produce Makefile.in
##
## # $Log$
## # Revision 1.3  2005/09/15 21:35:28  nadya
## # remove extra trailing slash
## #
## # Revision 1.2  2005/09/14 00:39:52  nadya
## # add mast sample sequences
## #
## # Revision 1.1.1.1  2005/07/25 23:21:55  nadya
## # Importing from meme-3.0.14, and adding configure/make
## #
##

AUTOMAKE_OPTIONS = no-dependencies

EXAMPLES = \
  logo_STRGGTCAN_+M00156.png \
  logo_STRGGTCAN_+M00515.png \
  logo_STRGGTCAN_+M01138.png \
  logo_STRGGTCAN_-M01282.png \
  logo_STRGGTCAN_+MA0159.1.png \
  logo_STRGGTCAN_+M00157.png \
  logo_STRGGTCAN_+M00528.png \
  logo_STRGGTCAN_+M01198.png \
  logo_STRGGTCAN_+MA0065.2.png \
  logo_STRGGTCAN_+MA0160.1.png \
  logo_STRGGTCAN_-M00172.png \
  logo_STRGGTCAN_+M00762.png \
  logo_STRGGTCAN_+M01202.png \
  logo_STRGGTCAN_+MA0066.1.png \
  logo_STRGGTCAN_+MA0258.1.png \
  logo_STRGGTCAN_-M00191.png \
  logo_STRGGTCAN_-M00763.png \
  logo_STRGGTCAN_-M01217.png \
  logo_STRGGTCAN_+MA0071.1.png \
  logo_STRGGTCAN_+MA0451.1.png \
  logo_STRGGTCAN_+M00202.png \
  logo_STRGGTCAN_+M00959.png \
  logo_STRGGTCAN_-M01227.png \
  logo_STRGGTCAN_+MA0072.1.png \
  logo_STRGGTCAN_-M00226.png \
  logo_STRGGTCAN_-M00965.png \
  logo_STRGGTCAN_+M01248.png \
  logo_STRGGTCAN_+MA0074.1.png \
  logo_STRGGTCAN_+M00242.png \
  logo_STRGGTCAN_+M01032.png \
  logo_STRGGTCAN_+M01268.png \
  logo_STRGGTCAN_-MA0112.1.png \
  logo_STRGGTCAN_+M00511.png \
  logo_STRGGTCAN_+M01033.png \
  logo_STRGGTCAN_+M01270.png \
  logo_STRGGTCAN_-MA0112.2.png \
  tomtom.xml \
  tomtom.html \
  tomtom.txt \
  STRGGTCAN.meme

if WEBSITE
exdir = $(WEBDIR)/doc/examples/tomtom_example_output_files
else
exdir = $(prefix)/doc/examples/tomtom_example_output_files
endif

ex_DATA = $(EXAMPLES)

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

EXTRA_DIST = $(EXAMPLES)
