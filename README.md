GimmeMotifs
===========



[![bioconda-badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io)
[![PyPI version](https://badge.fury.io/py/gimmemotifs.svg)](https://badge.fury.io/py/gimmemotifs)
[![Build Status](https://travis-ci.org/simonvh/gimmemotifs.svg?branch=master)](https://travis-ci.org/simonvh/gimmemotifs)
[![Code Health](https://landscape.io/github/simonvh/gimmemotifs/master/landscape.svg?style=flat)](https://landscape.io/github/simonvh/gimmemotifs/master)

Suite of motif tools, including a motif prediction pipeline for ChIP-seq experiments.

See [full GimmeMotifs documentation](http://gimmemotifs.readthedocs.org/) for detailed installation instructions and usage examples.

Quick start
-----------

Install using conda:

$ conda install gimmemotifs -c bioconda

Or with pip:

$ pip install gimmemotifs

Download a genome:

$ gimme genome hg38

Or index a genome directory with chromosome FASTA files on your computer:

$ gimme index /usr/share/genomes/hg19 hg19

Predict some motifs:

$ gimme motifs my_peaks.bed -g hg38 -n my_motifs 

