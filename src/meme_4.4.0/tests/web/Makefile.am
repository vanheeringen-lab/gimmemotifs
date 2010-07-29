## $Id: Makefile.am 3534 2008-12-03 20:24:51Z cegrant $
##
## Process this file with automake to produce Makefile.in
##

test-web: test-meme-cgi.pl
	./test-meme-cgi.pl $(URL)/cgi-bin/meme.cgi $(TEST_CONTACT)

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

EXTRA_DIST = $(DATAFILES)
