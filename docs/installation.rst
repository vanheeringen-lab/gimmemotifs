Installation
============

GimmeMotifs runs on Linux. On Windows 10 it will run fine using the `Windows Subsystem for Linux`_.
Mac OSX should work and is included in the build test. 
However, as I don't use it myself, unexpected issues might pop up. 
Let me know, so I can try to fix it.

.. _`Windows Subsystem for Linux`: https://docs.microsoft.com/en-us/windows/wsl/install-win10

.. _`Install GimmeMotifs`:

The easiest way to install
--------------------------

The preferred way to install GimmeMotifs is by using `conda
<https://docs.continuum.io/anaconda>`_. 
Activate the bioconda_ channel if you haven't used bioconda before.
You only have to do this once.

:: 

    $ conda config --add channels defaults
    $ conda config --add channels bioconda
    $ conda config --add channels conda-forge

You can install GimmeMotifs with one command. In the current environment:

::

    $ conda install gimmemotifs

Or create a specific environment:

::

    $ conda create -n gimme python=3 gimmemotifs
    
    # Activate the environment before you use GimmeMotifs
    $ conda activate gimme

GimmeMotifs only supports Python 3. Don't forget to activate the environment with ``source activate gimme`` whenever you want to use GimmeMotifs.

Installation successful? Good. Have a look at the :ref:`configuration<configuration>` section.

.. _`upgradegenome`:

Important note on upgrading from 0.11.1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The way genomes are installed and used has been changed from 0.11.1 to 0.12.0. 
Basically, we have switched to the faidx index used and supported by many other tools. 
This means that the old (<=0.11.1) GimmeMotifs index cannot be used by GimmeMotifs 0.12.0 and higher. 
You can re-install genomes using genomepy_, which is now the preferred tool for genome management for GimmeMotifs.
However, because of this change you can now also directly supply a genome FASTA instead of a genome name. 
Pre-indexing is not required anymore.

.. _bioconda: https://bioconda.github.io/
.. _genomepy: https://github.com/simonvh/genomepy

Alternative installation
------------------------

Prerequisites
+++++++++++++

These are the prerequisites for a full GimmeMotifs installation.

- bedtools http://bedtools.readthedocs.io
- UCSC genePredToBed http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToBed
- UCSC bigBedToBed http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed
- Perl + Algorithm::Cluster

In addition many of the motif tools (such as MEME) will need to be installed separately. Instructions for doing so are not included here.

Using pip
+++++++++

Installation from PyPI with ``pip`` is a relatively straightforward option.
Install with pip as follows:

:: 

    $ sudo pip install gimmemotifs

Or the (unstable) develop branch with the newest bells, whistles and bugs:

::

    $ sudo pip install git+https://github.com/vanheeringen-lab/gimmemotifs.git@develop

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

On first run GimmeMotifs will try to locate the tools you have
installed. If you have recently installed them, running an ``updatedb``
will be necessary. Using this option GimmeMotifs will create a
configuration file, the default is:

::

    ~/.config/gimmemotifs/gimmemotifs.cfg

This is a personal configuration file.

It is also possible to run the ``setup.py install`` command with the
``--prefix``, ``--home``, or ``--install-data`` options, to install in
GimmeMotifs in a different location (for instance, in your own home
directory). This should be fine, however, these alternative methods of
installing GimmeMotifs have not been extensively tested. 

.. _configuration:

Configuration
-------------

Genomes
+++++++

You will need genome FASTA files for a lot of the tools that are included 
with GimmeMotifs.

Download genomes automatically
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The most straightforward way to download and index a genome is to use
the ``genomepy`` tool, which is installed with GimmeMotifs.

::

    $ genomepy install hg19 UCSC --annotation

Here, the hg19 genome and accompanying gene annotation will be downloaded
from UCSC to the directory ``~/.local/share/genomes/hg19``. 
You can change this default location by creating/editing the file ``~/.config/genomepy/genomepy.yaml`` and change the following line:

::

    genome_dir: /data/genomes

Please note: in contrast to earlier versions of GimmeMotifs it is no longer necessary to index a genome.

Adding gene files
~~~~~~~~~~~~~~~~~

Note: If you used the ``genomepy`` command, annotation will be included automatically.

For some applications a gene file is used. This is a file containing gene
annotation in BED12 format. It should be located in the ``gene_dir``, 
which is defined in the configuration file (see below). 
The file needs to be named ``<index_name>.bed``, so for instance ``hg19.bed``.

.. _`other_configuration`:

Other configuration options
+++++++++++++++++++++++++++

All of GimmeMotifs configuration is stored in
``~/.config/gimmemotifs/gimmemotifs.cfg``. The
configuraton file is created at first run with  all defaults set,
but you can always edit it afterwards. It contains two sections ``main``
and ``params`` that take care of paths, file locations, parameter
settings etc. Additionally, every motif tool has it's own section. Let's
have a look at the options.

::

    [main]
    template_dir = /usr/share/gimmemotifs/templates
    score_dir = /usr/share/gimmemotifs/score_dists
    motif_databases = /usr/share/gimmemotifs/motif_databases
    gene_dir = /usr/share/gimmemotifs/genes
    tools_dir = /usr/share/gimmemotifs/tools

-  ``template_dir`` The location of the jinja2 html templates, used to
   generate the reports.

-  ``score_dir`` To generate p-values, a pre-calculated file with mean
   and sd of score distributions is needed. These are located here.

-  ``motif_databases`` For now contains only the JASPAR motifs.

-  ``gene_dir`` Directory with bed-files containing gene locations.
   This is needed to create promoter background sequences.

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
