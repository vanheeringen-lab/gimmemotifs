## $Id: Makefile.am 4911 2010-09-03 07:10:33Z james_johnson $
##
## Process this file with automake to produce Makefile.in
##
## # $Log$
## # Revision 1.3  2006/03/07 23:30:25  nadya
## # merge branches v3_5_1 and v3_5_2 back to the trunk
## #
## # Revision 1.2.6.1  2006/02/23 00:36:12  nadya
## # add gif images for styling with js and css
## #
## # Revision 1.2  2005/08/07 04:45:13  nadya
## # add missing host_logo.png
## #
## # Revision 1.1.1.1  2005/07/25 23:17:43  nadya
## # Importing from meme-3.0.14, and adding configure/make
## #
##

AUTOMAKE_OPTIONS = no-dependencies

MEME_IMAGES = \
  acb_logo.png \
  base.gif \
  bkg.jpg \
  blankpixel.gif \
	consensus-eq-1.png \
	consensus-eq-2.png \
  empty.gif \
	euclidean-eq.png \
  flowchart.gif \
  folder.gif \
  folderopen.gif \
	help.png \
  gomo.png \
  host_logo.png \
  join.gif \
  joinbottom.gif \
  line.gif \
  logo_alignment.png \
  mast.png \
  mcast-logo.fig \
  mcast-logo.png \
  meme.png \
  memesuite.png \
  memesuite.eps \
  memesuite.fig \
  memesuite_icon.png\
  memesuite_icon.ico\
  meta-meme.gif \
  minus.gif \
  minus.png \
  minusbottom.gif \
  page.gif \
  pasteur_logo.png \
	pearson-eq.png \
  plus.gif \
  plus.png \
  plusbottom.gif \
	sandelin-eq.png \
  sdsc.gif  

if WEBSITE
memeimgdir = $(WEBDIR)/images
memeimg_DATA = $(MEME_IMAGES)
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

EXTRA_DIST = $(MEME_IMAGES)
