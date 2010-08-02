# This script creates a RPM package of gimmemotifs for distribution
from subprocess import Popen
from shutil import copyfile

Popen(["python","setup.py","sdist"]).communicate()
copyfile("dist/gimmemotifs-0.50.tar.gz", "/root/rpmbuild/SOURCES/gimmemotifs-0.50.tar.gz")

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
