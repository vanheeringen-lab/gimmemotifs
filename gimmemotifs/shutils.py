# Copyright (c) 2009-2015 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

""" Odds and ends that for which I didn't (yet) find another place """

# Python imports
import os
import sys
from subprocess import Popen,PIPE

def run_command(cmd):
    #print args
    from subprocess import Popen
    p = Popen(cmd, shell=True)
    p.communicate()

def which(file):
    if not os.environ.has_key("PATH") or not os.environ["PATH"]:
        path = os.defpath
    else:
        path = os.environ["PATH"]

    for p in [os.path.join(x, file) for x in path.split(os.pathsep)]:
        if os.access(p, os.X_OK) and not os.path.isdir(p):
            return p

    p = Popen("locate %s" % file, shell=True, stdout=PIPE, stderr=PIPE)
    (stdout, stderr) = p.communicate()
    if not stderr:
        for p in stdout.split("\n"):
            if os.path.basename(p) == file and os.access(p, os.X_OK) and not os.path.isdir(p):
                return p

