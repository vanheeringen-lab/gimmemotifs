# This script creates a RPM package of gimmemotifs for distribution
import os
import sys
from subprocess import Popen, PIPE
from glob import glob
import shutil
from datetime import datetime

build_dir = "build/debian"
author = "Simon van Heeringen"
email = "simon.vanheeringen@gmail.com"
timestamp = datetime.now().strftime("%a, %d %b %Y %H:%M:%S +0100")

# Get package name
name = Popen(["python","setup.py","--name"], stdout=PIPE).communicate()[0].strip()
# Get package version
version  = Popen(["python","setup.py","--version"], stdout=PIPE).communicate()[0].strip()

package = "%s-%s" % (name, version)

Popen(["python","setup.py","sdist"]).communicate()
if not os.path.exists("build"):
	os.mkdir("build")

if os.path.exists("build/debian"):
	shutil.rmtree("build/debian")
os.mkdir("build/debian")
shutil.copyfile("dist/%s.tar.gz" % package, "%s/%s.tar.gz" % (build_dir, package))
Popen(["tar","xvzf","%s.tar.gz" % package], cwd=build_dir).communicate()
Popen(["dh_make -c bsd -e %s -s -f ../%s.tar.gz" % (email, package) ], cwd="%s/%s" % (build_dir, package), shell=True, stdin=PIPE).communicate("\n")

files = []
patterns = ["*ex", "*EX", "dirs", "docs", "README.Debian"]
for p in patterns:
	files += glob("%s/%s/debian/%s" % (build_dir, package, p))

for file in files:
	os.remove(file)

for file in ["rules"]:
	shutil.copyfile(file, ("%s/%s/debian/" % (build_dir, package)) + file)

if Popen(["python","--version"], stderr=PIPE).communicate()[1].find("2.7") != -1:
	shutil.copyfile("control2.7", "%s/%s/debian/control" % (build_dir, package))
	f = open("%s/%s/debian/pyversions" % (build_dir, package),"w")
	f.write("2.7-")
	f.close()

elif Popen(["python","--version"], stderr=PIPE).communicate()[1].find("2.6") != -1:
	shutil.copyfile("control2.6", "%s/%s/debian/control" % (build_dir, package))
	f = open("%s/%s/debian/pyversions" % (build_dir, package),"w")
	f.write("2.6-")
	f.close()
else:
	shutil.copyfile("control2.5", "%s/%s/debian/control" % (build_dir, package))
	f = open("%s/%s/debian/pyversions" % (build_dir, package),"w")
	f.write("2.5-")
	f.close()

f = open("%s/%s/debian/changelog" % (build_dir, package), "w")
f.write("gimmemotifs (%s) unstable; urgency=low\n" % version)
f.write("\n  * new version\n\n")
f.write(" -- %s <%s>  %s" % (author, email, timestamp)) 
f.close()

Popen(["debuild -us -uc"], cwd="%s/%s" % (build_dir, package), shell=True).communicate("")
