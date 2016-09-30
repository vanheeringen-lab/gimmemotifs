# GimmeMotifs

[![bioconda-badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io)
[![PyPI version](https://badge.fury.io/py/gimmemotifs.svg)](https://badge.fury.io/py/gimmemotifs)
[![Build Status](https://travis-ci.org/simonvh/gimmemotifs.svg?branch=master)](https://travis-ci.org/simonvh/gimmemotifs)
[![Code Health](https://landscape.io/github/simonvh/gimmemotifs/master/landscape.svg?style=flat)](https://landscape.io/github/simonvh/gimmemotifs/master)

Suite of motif tools, including a motif prediction pipeline for ChIP-seq experiments.

See [full GimmeMotifs documentation](http://gimmemotifs.readthedocs.org/) for detailed installation instructions and usage examples.

## Easy installation

The most straightforward way to install GimmeMotifs is via [conda](https://docs.continuum.io/anaconda/) using the [bioconda](https://bioconda.github.io/) channel.

`$ conda install gimmemotifs -c bioconda`

## Quick start

### Download a genome

Create a directory to store genome files.

`$ mkdir $HOME/genomes/`

To download and index a genome (all UCSC-supported genomes):

`$ gimme genome $HOME/genomes/ hg38`

Alternatively, you can index a genome directory with chromosome FASTA files on your computer.

`$ gimme index /usr/share/genomes/hg19 hg19`

### Predict some motifs:

`$ gimme motifs my_peaks.bed -g hg38 -n my_motifs`

## Frequently Asked Questions (FAQ)

#### I get the following error: "Invalid value for background argument".

Currently, this is a bug in the default configuration file. Run `gimme motifs` with the additional argument `-b gc,random`. 



