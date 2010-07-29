## $Id: Makefile.am 2599 2008-04-30 08:59:29Z tbailey $
##
## Process this file with automake to produce Makefile.in
##
## # $Log$
##

AUTOMAKE_OPTIONS = no-dependencies

EXAMPLES = \
  glam2scan.html 

if WEBSITE
exdir = $(WEBDIR)/doc/examples/glam2scan_example_output_files
else
exdir = $(prefix)/doc/examples/glam2scan_example_output_files
endif

ex_DATA = $(EXAMPLES)

# If .svnignore file doesn't exist
# create an empty one
.svnignore:
	touch .svnignore

# Update SVN svn:ignore property
# If .svnignore exists and longer then 0 bytes.
.svnignore.old:.svnignore
	if [ -s .svnignore ] ; \
  then \
    svn propset svn:ignore -F .svnignore .  ; \
    cp .svnignore .svnignore.old; \
  fi

EXTRA_DIST = $(EXAMPLES)
