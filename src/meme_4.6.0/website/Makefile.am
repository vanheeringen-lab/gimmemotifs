## $Id: Makefile.am 5296 2011-01-11 03:53:36Z tbailey $
##
## Process this file with automake to produce Makefile.in
##
#

AUTOMAKE_OPTIONS = no-dependencies

if WEBSITE
SUBDIRS = html cgi-bin images scripts

## If we are building the web application, build each
## of the child subdirectories.
# Make some directories needed by the web application
install-data-local:
	mkdir -p $(WEBDIR)/output
	chmod a+w $(WEBDIR)/output
	touch $(WEBDIR)/output/index.html	# prevent dir listing
	mkdir -p $(WEBDIR)/logs
	chmod a+w $(WEBDIR)/logs
	mkdir -p $(WEBDIR)/queue
	chmod a+w $(WEBDIR)/queue
	mkdir -p $(WEBDIR)/user-files
	chmod a+w $(WEBDIR)/user-files

else
if WEBSERVICE
SUBDIRS = scripts cgi-bin
else
SUBDIRS = scripts 
endif
endif

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


