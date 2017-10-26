Installation
============

GimmeMotifs runs on Linux. Definitely not on Windows, sorry. 
Mac OSX should work and is included in the build test. 
However, as I don't use it myself, unexpected issues might pop up. 
Let me know, so I can try to fix it.

.. _`Install GimmeMotifs`:

The easiest way to install
--------------------------

The preferred way to install GimmeMotifs is by using `conda
<https://docs.continuum.io/anaconda>`_. 
Activate the bioconda_ channel if you haven't used bioconda before.
You only have to do this once.

:: 

    $ conda config --add channels r
    $ conda config --add channels defaults
    $ conda config --add channels conda-forge
    $ conda config --add channels bioconda

You can install GimmeMotifs with one command. In the current environment:

::

    $ conda install gimmemotifs

Or create a specific environment:

::

    $ conda create -n gimme python gimmemotifs
    
    # Activate the environment before you use GimmeMotifs
    $ source activate gimme

Python 3 is the preferred version, however, GimmeMotifs also supports Python 2.
Don't forget to activate the environment with ``source activate gimme`` whenever
you want to use GimmeMotifs.

Installation successful? Good. Have a look at the :ref:`configuration<configuration>` section.

.. _bioconda: https://bioconda.github.io/

Alternative installation
------------------------

Prerequisites
+++++++++++++

These are the prerequisites for a full GimmeMotifs installation.

- bedtools http://bedtools.readthedocs.io
- UCSC genePredToBed http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToBed
- UCSC bigBedToBed http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed
- R + RobustRankAggreg https://cran.r-project.org/web/packages/RobustRankAggreg/index.html
- Perl + Algorithm::Cluster

Using pip
+++++++++

Installation from PyPI with ``pip`` is a relatively straightforward option.
Install with pip as follows:

:: 

    $ sudo pip install gimmemotifs

Or the (unstable) develops branch with the newest bells, whistles and bugs:

::

    $ sudo pip install git+https://github.com/simonvh/gimmemotifs.git@develop

If you don't have root access, see the option below.

Using pip in a virtualenv
+++++++++++++++++++++++++

Ubuntu prerequisites
~~~~~~~~~~~~~~~~~~~~

To install GimmeMotifs in a virtualenv, several Python packages need to be built from source. 

Install the necessary packages to build numpy, scipy, matplotlib and GimmeMotifs:

::

    sudo apt-get install python-pip python-dev build-essential libatlas-base-dev \
    gfortran liblapack-dev libatlas-base-dev cython libpng12-dev libfreetype6-dev \
    libgsl0-dev

Install via pip
~~~~~~~~~~~~~~~

Create a virtualenv and activate it according to the 
`documentation
<https://virtualenv.readthedocs.org/en/latest/userguide.html#usage>`_.

Install numpy:

::

    $ pip install numpy


Now you can install GimmeMotifs using pip. Latest stable release:

::

    $ pip install gimmemotifs


Installation from source
++++++++++++++++++++++++

Did I mention conda? 

You know bioconda is amazing, right?

So...


These instructions are not up-to-date! Basically, you're on your own!

Make sure to install all required dependencies.

You can download the lastest stable version of GimmeMotifs at:

| https://github.com/simonvh/gimmemotifs/releases

Start by unpacking the source archive

::

    tar xvzf gimmemotifs-0.11.0.tar.gz
    cd gimmemotifs-0.11.0

You can build GimmeMotifs with the following command:

::

    python setup.py build

Run the tests to check if the basics work correctly:

::

    python run_tests.py

If you encounter no errors, go ahead with installing GimmeMotifs (root
privileges required):

::

    sudo python setup.py install

During installation GimmeMotifs will try to locate the tools you have
installed. If you have recently installed them, running an ``updatedb``
will be necessary. Using this option GimmeMotifs will create a
configuration file, the default is:

::

    /usr/share/gimmemotifs/gimmemotifs.cfg

This is a system-wide configuration that can be used by all users.

It is also possible to run the ``setup.py install`` command with the
``–prefix``, ``–home``, or ``–install-data`` options, to install in
GimmeMotifs in a different location (for instance, in your own home
directory). This should be fine, however, these alternative methods of
installing GimmeMotifs have not been extensively tested. Please note
that in this case the configuration file will be created, but every user
will have to put this configuration file in his/her home directory:
``~/.gimmemotifs.cfg``. The install script will also inform you of this
during install.  


.. _configuration:

Configuration
-------------

Genomes
+++++++

You will need genome FASTA files for a lot of the tools that are included 
with GimmeMotifs.

Download from UCSC
~~~~~~~~~~~~~~~~~~

The most straightforward way to download and index a genome is to use
the ``gimme genome`` tool.

::

    $ gimme genome $HOME/genomes hg19

Here, the hg19 genome and accompanying gene annotation will be downloaded
from UCSC to the directory ``$HOME/genomes/hg19``. 
This should work for all genomes supported by UCSC. 

Index a genome
~~~~~~~~~~~~~~

Alternatively, you can index a set of genome FASTA files that you already
have locally. The FASTA files should be organized in one
directory with *one file per chromosome or scaffold*, with the filename
being the chromosome name with an extension of ``.fa``, ``.fsa`` or
``.fasta``. Then you can run the following command:

::

    gimme index /dir/to/fasta/files/ genome_name

For instance, if I wanted to index the human genome (version hg19) on my
computer, where all fasta files are located in the directory
``/usr/share/genome/hg19`` I would run the following command:

::

    gimme index /usr/share/genome/hg19/ hg19

**Note: if you installed GimmeMotifs as root, the** ``gimme index`` **command
will need to be run as root too** 

Adding gene files
~~~~~~~~~~~~~~~~~

For some applications a gene file is used. This is a file containing gene
annotation in BED12 format. It should be located in the ``gene_dir``, 
which is defined in the configuration file (see below). 
The file needs to be named ``<index_name>.bed``, so for instance ``hg19.bed``.
If you used the ``gimme genome`` command, 
annotation will be included automatically.

.. _`other_configuration`:

Other configuration options
+++++++++++++++++++++++++++

All of GimmeMotifs configuration is stored in
``/usr/share/gimmemotifs/gimmemotifs.cfg`` or ``~/.gimmemotifs.cfg``. If
the file ``~/.gimmemotifs.cfg`` exists in your home directory this will
always have precedence over the system-wide configuration. The
configuraton file is created at installation time with all defaults set,
but you can always edit it afterwards. It contains two sections ``main``
and ``params`` that take care of paths, file locations, parameter
settings etc. Additionally, every motif tool has it’s own section. Let’s
have a look at the options.

::

    [main]
    index_dir = /usr/share/gimmemotifs/genome_index
    template_dir = /usr/share/gimmemotifs/templates
    seqlogo = /usr/local/bin/seqlogo
    score_dir = /usr/share/gimmemotifs/score_dists
    motif_databases = /usr/share/gimmemotifs/motif_databases
    gene_dir = /usr/share/gimmemotifs/genes
    tools_dir = /usr/share/gimmemotifs/tools

-  ``index_dir`` The location of the indeces of the genome fasta-files.

-  ``template_dir`` The location of the KID html templates, used to
   generate the reports.

-  ``seqlogo`` The seqlogo executable.

-  ``score_dir`` To generate p-values, a pre-calculated file with mean
   and sd of score distributions is needed. These are located here.

-  ``motif_databases`` For now contains only the JASPAR motifs.

-  ``gene_dir`` Directory with bed-files containing gene locations for
   every indexed organism. This is needed to create the matched genomic
   background.

-  ``tools_dir`` Here all tools included with GimmeMotifs are stored.

::

    [params]
    fraction = 0.2
    use_strand = False
    abs_max = 1000
    analysis = medium
    enrichment = 1.5
    width = 200
    lwidth = 500
    genome = hg19
    background = gc,random
    cluster_threshold = 0.95
    available_tools = MDmodule,MEME,Weeder,GADEM,MotifSampler,trawler,Improbizer,BioProspector,Posmo,ChIPMunk,JASPAR,AMD,HMS,Homer
    tools = MDmodule,MEME,Weeder,MotifSampler,trawler,Improbizer,BioProspector,Posmo,ChIPMunk,JASPAR,AMD,HMS,Homer
    pvalue = 0.001
    max_time = None
    ncpus = 2
    motif_db = gimme.vertebrate.v3.1.pwm
    scan_cutoff = 0.9
    use_cache = False
    markov_model = 1
    
This section specifies all the default GimmeMotifs parameters. Most of
these can also be specified at the command-line when running
GimmeMotifs, in which case they will override the parameters specified

Configuration of MotifSampler
+++++++++++++++++++++++++++++

If you want to use MotifSampler there is one more step that you'll have
to take *after* installation of GimmeMotifs. For every organism, you will
need a MotifSampler background. Note that human (hg19, hg38) and mouse (mm9, mm10) background models are included, so for these
organisms MotifSampler will work out of the box. For other organisms the
necessary background files can be created with ``CreateBackgroundModel``
(which is included with GimmeMotifs or can be downloaded from the same
site as MotifSampler). The background model file needs to be saved in
the directory ``/usr/share/gimmemotifs/MotifSampler`` and it should be
named ``<organism_index_name>.bg``. So, for instance, if I downloaded
the human epd background
(``epd_homo_sapiens_499_chromgenes_non_split_3.bg``), this file should
be saved as ``/usr/share/gimmemotifs/MotifSampler/hg19.bg``.
here.
