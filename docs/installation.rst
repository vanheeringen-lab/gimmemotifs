Installation
============

GimmeMotifs runs on Linux. On Windows it will run fine using the `Windows Subsystem for Linux`_.

..  NOTE: nope. it hasn't worked in a while.
    Mac OSX should work and is included in the build test.
    However, as I don't use it myself, unexpected issues might pop up.
    Let me know, so I can try to fix it.

.. _`Windows Subsystem for Linux`: https://docs.microsoft.com/en-us/windows/wsl/install-win10

.. _`Install GimmeMotifs`:



Conda - the easy way
--------------------

The preferred way to install GimmeMotifs is by using conda_.

.. _conda: https://docs.continuum.io/anaconda
.. _mamba: https://github.com/mamba-org/mamba
.. _bioconda: https://bioconda.github.io/
.. _271: https://github.com/vanheeringen-lab/gimmemotifs/issues/271

:: 

    $ conda config --add channels bioconda
    $ conda config --add channels conda-forge

You can install GimmeMotifs with one command. In the current environment:

::

    $ conda install gimmemotifs

Or create a specific environment:

::

    $ conda create -n gimme gimmemotifs
    
    # Activate the environment before you use GimmeMotifs
    $ conda activate gimme

Installation successful? Good. Have a look at the :ref:`configuration<configuration>` section.



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

    # installs gimmemotifs (in editable mode)
    $ pip install --no-deps --no-cache-dir --use-pep517 -v -e .

    # test if the install was successful
    $ gimme -h

Once installed, you can edit the code in the `gimmemotifs` folder, and the changes are immediately active!
Check out how good your fixes are with unit tests:

::

    $ pytest -vvv --disable-pytest-warnings



Configuration
=============

.. _configuration:
.. _`other_configuration`:

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



MotifSampler
------------

If you want to use MotifSampler there is one more step that you'll have to take *after* installation of GimmeMotifs.
For every organism, you will need a MotifSampler background.
Note that human (hg19, hg38) and mouse (mm9, mm10) background models are included, so for these organisms MotifSampler will work out of the box.
For other organisms the necessary background files can be created with ``CreateBackgroundModel`` (which is included with GimmeMotifs or can be downloaded from the same site as MotifSampler).
The background model file needs to be saved in the directory ``~/.share/gimmemotifs/MotifSampler`` and it should be named ``<organism_index_name>.bg``.
So, for instance, if I downloaded the human epd background (``epd_homo_sapiens_499_chromgenes_non_split_3.bg``), this file should be saved as ``~/.share/gimmemotifs/MotifSampler/hg19.bg`` here.
