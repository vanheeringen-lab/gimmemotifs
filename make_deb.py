# This script creates a RPM package of gimmemotifs for distribution
import os
from subprocess import Popen, PIPE
from glob import glob
import shutil

Popen(["python","setup.py","sdist"]).communicate()
if not os.path.exists("build"):
	os.mkdir("build")

if os.path.exists("build/debian"):
	shutil.rmtree("build/debian")
os.mkdir("build/debian")
shutil.copyfile("dist/gimmemotifs-0.50.tar.gz", "build/debian/gimmemotifs-0.50.tar.gz")
Popen(["tar","xvzf","gimmemotifs-0.50.tar.gz"], cwd="build/debian").communicate()
Popen(["dh_make -c bsd -s -f ../gimmemotifs-0.50.tar.gz"], cwd="build/debian/gimmemotifs-0.50", shell=True, stdin=PIPE).communicate("\n")

files = glob("build/debian/gimmemotifs-0.50/debian/*.ex") + glob("build/debian/gimmemotifs-0.50/debian/*.EX") + ["build/debian/gimmemotifs-0.50/debian/dirs", "build/debian/gimmemotifs-0.50/debian/docs", "build/debian/gimmemotifs-0.50/debian/README.Debian"]

for file in files:
	os.remove(file)

for file in ["control", "pyversions", "rules"]:
	shutil.copyfile(file, "build/debian/gimmemotifs-0.50/debian/" + file)	
Popen(["debuild -us -uc"], cwd="build/debian/gimmemotifs-0.50", shell=True).communicate("")

#Popen(["dh_make]).communicate()


#lines = [line for line in open(SPECFILE)]
#f = open(SPECFILE, "w")
#f.write("%define debug_package %{nil}\n")
#f.write("%define __spec_install_port /usr/lib/rpm/brp-compress\n")
#for line in lines:
#	f.write(line)
#f.close()

#Popen(["rpmbuild","-ba", "dist/gimmemotifs.spec"])
