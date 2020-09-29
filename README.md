# GimmeMotifs

[![bioconda-badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io)
[![PyPI version](https://badge.fury.io/py/gimmemotifs.svg)](https://badge.fury.io/py/gimmemotifs)
[![Build Status](https://travis-ci.org/vanheeringen-lab/gimmemotifs.svg?branch=master)](https://travis-ci.org/vanheeringen-lab/gimmemotifs)
[![Documentation Status](https://readthedocs.org/projects/gimmemotifs/badge/?version=master)](http://gimmemotifs.readthedocs.io/en/master/?badge=master)
[![Maintainability](https://api.codeclimate.com/v1/badges/43b1243e35b8a00fac8a/maintainability)](https://codeclimate.com/github/vanheeringen-lab/gimmemotifs/maintainability)
[![Test Coverage](https://api.codeclimate.com/v1/badges/43b1243e35b8a00fac8a/test_coverage)](https://codeclimate.com/github/vanheeringen-lab/gimmemotifs/test_coverage)

[![Anaconda-Server Badge](https://anaconda.org/bioconda/gimmemotifs/badges/downloads.svg)](https://anaconda.org/bioconda/gimmemotifs)
[![DOI](https://zenodo.org/badge/676678.svg)](https://zenodo.org/badge/latestdoi/676678)

Suite of motif tools, including a motif prediction pipeline for ChIP-seq experiments.

See [full GimmeMotifs documentation](http://gimmemotifs.readthedocs.org/) for detailed installation instructions and usage examples.

For documentation on the development version see [here](http://gimmemotifs.readthedocs.org/en/latest/).

The manuscript describing this latest release is available on [biorRxiv](https://doi.org/10.1101/474403) as a preprint and can be cited as:

> [**GimmeMotifs: an analysis framework for transcription factor motif analysis**](https://doi.org/10.1101/474403) <br>
Niklas Bruse, Simon J. van Heeringen<br>
_bioRxiv_ (2018) DOI: [10.1101/474403](https://doi.org/10.1101/474403)

You can interactively try out the Python API in a Jupyter notebook using binder: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/vanheeringen-lab/gimmemotifs/develop?filepath=docs%2Fapi_examples.ipynb)

## We need your help!

GimmeMotifs was originally developed for our own needs but we would really like it to be useful to the wider community. However, this also depends on your input. Let us know what you think! What features are missing? Which tutorial would you like to see? What part of the documentation is unclear? Have great ideas for future developments? Maybe you even want to join in developing this software?

[Let us know!](https://github.com/simonvh/gimmemotifs/issues/)



## Easy installation

The most straightforward way to install GimmeMotifs is via [conda](https://docs.continuum.io/anaconda/) using the [bioconda](https://bioconda.github.io/) channel.

If you have not used bioconda before, first set up the necessary channels (in this order!). You only have to do this once.

```
$ conda config --add channels defaults
$ conda config --add channels bioconda
$ conda config --add channels conda-forge
```

You can now install GimmeMotifs with one command:

```
# Create an environment called gimme with all dependencies
$ conda create -n gimme python=3 gimmemotifs

# Activate the environment
$ conda activate gimme
```

Python 3 is the required, from version 0.13.0 on GimmeMotifs no longer supports Python 2. 
Don't forget to activate the environment with `conda activate gimme` whenever you want to use GimmeMotifs.

## Quick start

### Predict some de novo motifs:

```
$ gimme motifs my_peaks.bed my_motifs -g /data/genomes/hg38/hg38.fa --denovo
```

### Download a genome

The example above assumes that you have the hg38 genome in
`/data/genomes/hg38/hg38.fa`. 
GimmeMotifs can also use genomes installed by
[genomepy](http://github.com/simonvh/genomepy).

You can configure the directory where genomepy stores genomes by editing
`~/.config/genomepy/genomepy.yaml`

``` 
genome_dir: /data/genomes
``` 

To download a genome from UCSC:

```
$ genomepy install hg38 --annotation  # genomepy >=0.9.0
```

Now you can specify this genome for GimmeMotifs by name.

```
$ gimme motifs my_peaks.bed -g hg38 -n my_motifs
```

## Help 


* Full documentation:
  [http://gimmemotifs.readthedocs.io/](http://gimmemotifs.readthedocs.io/).
* Check the [FAQ](http://gimmemotifs.readthedocs.io/en/master/faq.html#faq) for
  common issues.
* The preferred way to get support is through the Github
  [issues](https://github.com/simonvh/gimmemotifs/issues/) page
* Finally, you can reach me by [mail](simon.vanheeringen@gmail.com) or
  [Twitter](https://twitter.com/svheeringen).


