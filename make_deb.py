# This script creates a RPM package of gimmemotifs for distribution
import os
import sys
from subprocess import Popen, PIPE
from glob import glob
import shutil

build_dir = "build/debian"
email = "s.vanheeringen@ncmls.ru.nl"

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

for file in ["control", "pyversions", "rules"]:
	shutil.copyfile(file, ("%s/%s/debian/" % (build_dir, package)) + file)	
Popen(["debuild -us -uc"], cwd="%s/%s" % (build_dir, package), shell=True).communicate("")
