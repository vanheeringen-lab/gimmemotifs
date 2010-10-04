# This script creates a RPM package of gimmemotifs for distribution
from subprocess import Popen,PIPE
from shutil import copyfile
import os

# Get package name
name = Popen(["python","setup.py","--name"], stdout=PIPE).communicate()[0].strip()
# Get package version
version  = Popen(["python","setup.py","--version"], stdout=PIPE).communicate()[0].strip()

package = "%s-%s" % (name, version)

Popen(["python","setup.py","sdist"]).communicate()
copyfile("dist/%s.tar.gz" % (package), os.path.expanduser("~/rpmbuild/SOURCES/%s.tar.gz" % package))

SPECFILE="dist/gimmemotifs.spec"
Popen(["python","setup.py","bdist_rpm","--spec-only"]).communicate()
lines = [line for line in open(SPECFILE)]
f = open(SPECFILE, "w")
f.write("%define debug_package %{nil}\n")
f.write("%define __spec_install_port /usr/lib/rpm/brp-compress\n")
for line in lines:
	f.write(line)
f.close()

Popen(["rpmbuild","-ba", "dist/gimmemotifs.spec"])
