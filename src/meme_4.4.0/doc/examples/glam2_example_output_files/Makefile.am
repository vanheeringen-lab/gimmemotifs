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
  glam2.html \
  logo1.png \
  logo_ssc1.png \
  logo2.png \
  logo_ssc2.png \
  logo3.png \
  logo_ssc3.png \
  logo4.png \
  logo_ssc4.png \
  logo5.png \
  logo_ssc5.png \
  logo6.png \
  logo_ssc6.png \
  logo7.png \
  logo_ssc7.png \
  logo8.png \
  logo_ssc8.png \
  logo9.png \
  logo_ssc9.png 

if WEBSITE
exdir = $(WEBDIR)/doc/examples/glam2_example_output_files
else
exdir = $(prefix)/doc/examples/glam2_example_output_files
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
