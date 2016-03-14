#!/usr/bin/env python
# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

from gimmemotifs.genome_index import *
from gimmemotifs.config import *
from gimmemotifs.utils import is_valid_bedfile
from optparse import OptionParser
import sys
import os
import shutil

parser = OptionParser()
parser.add_option("-n", "--name", dest="name", help="Organism name to use (for example: hg18)", metavar="NAME")
parser.add_option("-f", "--fastadir", dest="fasta_dir", help="Directory containing fastafiles", metavar="DIR")
parser.add_option("-g", "--genes", dest="genes", help="Genes in BED format", metavar="FILE")
parser.add_option("-r", "--replace_index", dest="replace", help="Replac/overwrite index if it already exists", default=False, action="store_true")
	
(options, args) = parser.parse_args()

if not options.fasta_dir or not options.name or not options.genes:
	parser.print_help()
	sys.exit(1)
	
config = MotifConfig()
fasta_dir = options.fasta_dir
index_dir = config.get_index_dir()
name = options.name
gene_file = options.genes

# Check validity gene file
if not is_valid_bedfile(gene_file):
	sys.stderr.write("Gene file is not a valid BED file")
	sys.exit()

# Check if index already exist
index_dir = os.path.join(index_dir, options.name)
if os.path.exists(index_dir) and not options.replace:
	print "Index with name %s already exists. Specify -r option to replace this index.\n" % name
	sys.exit()

# Copy gene file to gene directory
gene_dir = config.get_gene_dir()
try:
	shutil.copy(gene_file, gene_dir)
except (shutil.Error,IOError), e:
	if e.args[0] == 13:
		print "No permission to copy gene file in directory %s. Superuser access needed?" % gene_dir
		sys.exit()
	elif e.args[0].endswith("same file"):
		pass
	else:
		print e
		sys.exit()

# Create the index
g = GenomeIndex()
g = g.create_index(fasta_dir, index_dir)

print
print "Genome %s indexed and ready for use with GimmeMotifs" % name
print
if config.is_configured("MotifSampler"):
	print "To use MotifSampler, an additional step is necessary. Please provide a MotifSampler background file called %s.MotifSampler.bg in the directory %s" % (name, config.get_bg_dir())
