prefix = /usr/local
exec_prefix = ${prefix}
bindir = ${exec_prefix}/bin

all:    run.csh

run.csh:
	echo "#!/bin/csh" > $@
	echo >> $@
	echo "$(bindir)/gadem -fseq ./p53KnownSitesEmbeddedInSimulated.seq -verbose 1" >> $@
	chmod a+x $@

distclean:
	rm -f run.csh

clean:
	rm -f run.csh

install:
