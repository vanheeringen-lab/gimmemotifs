Installation
============

GimmeMotifs runs on Linux. On Windows 10 it will run fine using the `Windows Subsystem for Linux`_.

..  NOTE: nope. it hasn't worked in a while.
    Mac OSX should work and is included in the build test.
    However, as I don't use it myself, unexpected issues might pop up.
    Let me know, so I can try to fix it.

.. _`Windows Subsystem for Linux`: https://docs.microsoft.com/en-us/windows/wsl/install-win10

.. _`Install GimmeMotifs`:

Conda - the easy way
--------------------

The preferred way to install GimmeMotifs is by using conda_.
Activate the required channels and install mamba_ (you only have to do this once).

In this example, conda and mamba versions are pinned due to a bug with mamba.
For more information, see issue 271_.

.. _conda: https://docs.continuum.io/anaconda
.. _mamba: https://github.com/mamba-org/mamba
.. _bioconda: https://bioconda.github.io/
.. _271: https://github.com/vanheeringen-lab/gimmemotifs/issues/271

:: 

    $ conda config --add channels defaults
    $ conda config --add channels bioconda
    $ conda config --add channels conda-forge
    $ conda install -c conda-forge "conda>=4.12" "mamba>=0.27"

You can install GimmeMotifs with one command. In the current environment:

::

    $ mamba install gimmemotifs

Or create a specific environment:

::

    $ mamba create -n gimme gimmemotifs
    
    # Activate the environment before you use GimmeMotifs
    $ mamba activate gimme

Installation successful? Good. Have a look at the :ref:`configuration<configuration>` section.

.. _upgradegenome:

Upgrading from 0.11.1
^^^^^^^^^^^^^^^^^^^^^

The way genomes are installed and used has been changed from 0.11.1 to 0.12.0.
Basically, we have switched to the faidx index used and supported by many other tools.
This means that the old (<=0.11.1) GimmeMotifs index cannot be used by GimmeMotifs 0.12.0 and higher.
You can re-install genomes using genomepy_, which is now the preferred tool for genome management for GimmeMotifs.
However, because of this change you can now also directly supply a genome FASTA instead of a genome name.
Pre-indexing is not required anymore.

.. _genomepy: https://github.com/vanheeringen-lab/genomepy

..  NOTE: abbreviated
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

Pip
---

Installation from PyPI with ``pip`` is a relatively straightforward option.
Install with pip as follows:

::

    $ pip install gimmemotifs

Or the (unstable) develop branch with the newest bells, whistles and bugs:

::

    $ pip install git+https://github.com/vanheeringen-lab/gimmemotifs.git@develop

Note that several dependencies and many of the motif tools (such as MEME) need to be installed separately.
Instructions for doing so are not included here.

..  NOTE: Lets keep it simple, with Conda, PIP or Source
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


Source - developers install
---------------------------

Want to fix that darned bug yourself?
Want to try out the latest features?

Well look no further!
You can install the develop branch with the newest bells, whistles and bugs:

::

    # download the gimmemotifs code
    $ git clone https://github.com/vanheeringen-lab/gimmemotifs.git
    $ cd gimmemotifs
    $ git checkout develop

    # setup the gimme conda environment
    $ conda env create -f requirements.yaml
    $ conda activate gimme
    $ python setup.py build  # installs the motif discovery tools
    $ pip install -e .       # installs gimmemotifs (in editable mode)

    # test if the install was successful
    $ gimme -h

Once installed, you can edit the code in the `gimmemotifs` folder, and the changes are immediately active!
Check out how good your fixes are with unit tests:

::

    $ pytest -vvv --disable-pytest-warnings

.. NOTE: I've replaced this with the editable install
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
=============

.. _`other_configuration`:

The configuration file
----------------------

All of GimmeMotifs' configuration is stored in ``~/.config/gimmemotifs/gimmemotifs.cfg``.
The configuration file is created at first run with all defaults set, but you can always edit it afterwards.
It contains two sections ``main`` and ``params`` that take care of paths, file locations, parameter settings etc.
Additionally, every motif tool has it's own section.
Let's have a look at the options.

::

    [main]
    bg = bg
    template_dir = templates
    score_dir = score_dists
    gene_dir = genes
    motif_databases = motif_databases
    tools = included_tools/

-  ``template_dir`` The location of the jinja2 html templates, used to
   generate the reports.

-  ``score_dir`` To generate p-values, a pre-calculated file with mean
   and sd of score distributions is needed. These are located here.

-  ``gene_dir`` Directory with bed-files containing gene locations.
   This is needed to create promoter background sequences.

-  ``motif_databases`` Contains various motif databases.

-  ``tools`` Here all tools included with GimmeMotifs are stored.

::

    [params]
    fraction = 0.2
    use_strand = False
    abs_max = 1000
    analysis = xl
    enrichment = 1.5
    size = 200
    lsize = 500
    background = gc,random
    cluster_threshold = 0.95
    scan_cutoff = 0.9
    available_tools = AMD,BioProspector,ChIPMunk,DiNAMO,GADEM,HMS,Homer,Improbizer,MDmodule,MEME,MEMEW,MotifSampler,Posmo,ProSampler,Trawler,Weeder,XXmotif,Yamda
    tools = BioProspector,Homer,MEME
    pvalue = 0.001
    max_time = -1
    ncpus = 12
    motif_db = gimme.vertebrate.v5.0.pfm
    use_cache = False

This section specifies all the default GimmeMotifs parameters. Most of
these can also be specified at the command-line when running
GimmeMotifs, in which case they will override the parameters specified.

Input Data
==========

Genomes - and how to get them
-----------------------------

You will need genome FASTA files for a lot of the tools that are included with GimmeMotifs.

The most straightforward way to download and index a genome is to use the ``genomepy`` tool, which is installed with GimmeMotifs.

::

    $ genomepy install hg38 --provider UCSC --annotation

Here, the hg38 genome and accompanying gene annotation will be downloaded from UCSC to the directory ``~/.local/share/genomes/hg38``.
You can change this default location by editing the file ``~/.config/genomepy/genomepy.yaml`` and change the following line:

::

    genomes_dir: /data/genomes

If this file does not exist, you can generate it with ``genomepy config generate``.
After downloading a genome with genomepy, you can use its name (e.g. ``hg38``) for gimme commands.

.. I think this is outdated:
    Adding gene annotation files
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    Note: If you used the ``genomepy`` command, annotation will be included automatically.

    For some applications a gene file is used. This is a file containing gene
    annotation in BED12 format. It should be located in the ``gene_dir``,
    which is defined in the configuration file (see below).
    The file needs to be named ``<index_name>.bed``, so for instance ``hg19.bed``.

MotifSampler
------------

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
