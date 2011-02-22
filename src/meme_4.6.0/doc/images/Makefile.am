AUTOMAKE_OPTIONS = no-dependencies

%.png: %.fig
	fig2dev -L png -S 4 $< > $@

MEME_IMAGES = \
  cbrc_logo.png \
  fimo_logo.png \
  glam2_logo.png \
  glam2scan_logo.png \
  mcast-logo.png \
  memechip_logo.png \
  nbcr_logo_trans.gif \
  ncrr_logo.gif \
  schematic.fig \
  schematic.png \
  sdsc_logo.png \
  spamo_logo.png \
  tomtom_logo.png \
  uq_logo.png \
  ucsd_logo.png \
  UWlogo150p.gif \
  webserver_schema.png


if WEBSITE
memeimgdir = $(WEBDIR)/doc/images
else
memeimgdir = $(prefix)/doc/images
endif

memeimg_DATA = $(MEME_IMAGES)

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

EXTRA_DIST = $(MEME_IMAGES)
