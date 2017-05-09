# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.
""" Odds and ends that for which I didn't (yet) find another place."""
# Python imports
import os
import subprocess as sp

def which(fname):
    """Find location of executable."""
    if "PATH" not in os.environ or not os.environ["PATH"]:
        path = os.defpath
    else:
        path = os.environ["PATH"]

    for p in [os.path.join(x, fname) for x in path.split(os.pathsep)]:
        if os.access(p, os.X_OK) and not os.path.isdir(p):
            return p

    p = sp.Popen("locate %s" % fname, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    (stdout, stderr) = p.communicate()
    if not stderr:
        for p in stdout.decode().split("\n"):
            if (os.path.basename(p) == fname) and (
                os.access(p, os.X_OK)) and (
                not os.path.isdir(p)):
                return p

def find_by_ext(dirname, ext):
    """Find all files in a directory by extension."""
    # Get all fasta-files
    try: 
        files = os.listdir(dirname) 
    except OSError: 
        if os.path.exists(dirname): 
            cmd = "find {0} -maxdepth 1 -name \"*\"".format(dirname) 
            p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE) 
            stdout, _stderr = p.communicate() 
            files = [os.path.basename(fname) for fname in stdout.decode().splitlines()] 
        else: 
            raise 
     
    retfiles = [os.path.join(dirname, fname) for fname in files if 
                    os.path.splitext(fname)[-1] in ext] 
 
    return retfiles 
