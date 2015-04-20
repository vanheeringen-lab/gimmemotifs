Configuration
=============

Data sources
------------

You will need some genome fasta files for any motif-prediction if you
want to run GimmeMotifs with BED files as input (which is recommended,
although GimmeMotifs will also work on FASTA files). To get from a BED
file to the sequence information these genomic fasta files are
absolutely required. The fasta files should be organized in one
directory with one file per chromosome or scaffold, with the filename
being the chromosome name with an extension of ``.fa``, ``.fsa`` or
``.fasta``. No exceptions, no different layouts. A good source is the
UCSC Genome Browser database . For instance, the human hg19 files needed
to run the examples included with GimmeMotifs can be downloaded here:

::

    ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.zip
     

All fasta files need to be indexed before GimmeMotifs can use them, see
section :ref:`indexing`.


.. _indexing:

Indexing the genomes
--------------------

All the genomes that you want to use with GimmeMotifs will need to be
indexed for (relatively) fast retrieval of sequences. You can do this,
once you have installed GimmeMotifs, by running the following command
(as root or with sudo):

::

    gimme index /dir/to/fasta/files/ genome_name

For instance, if I wanted to index the human genome (version hg19) on my
computer, where all fasta files are located in the directory
``/usr/share/genome/`` I would run the following command:

::

    gimme index /usr/share/genome/hg19/ hg19

Repeat this step for every additional genome or organism that you want
to use GimmeMotifs with. Please note: for Weeder, currently only hg18,
hg19, mm9, rn4, dm3, fr2, danRer6, danRer7, galGal3, ce3, anoGam1,
yeast, sacCer2, xenTro2 and xenTro3 are supported as organism names
(following the UCSC naming convention). This will be fixed as a
configuration file in a later release.

Adding gene files
-----------------

When using the ``genomic_matched`` background setting (which is the
default), there needs to be a file describing genes in BED format in the
``gene_dir``, which is defined in the configuration file. By default
this is: ``/usr/share/gimmemotifs/genes/``. The file needs to be named
``<index_name>.bed``, so for instance ``hg19.bed``. By default
``hg18.bed``, ``hg19.bed``, ``mm9.bed`` and ``xenTro2.bed`` are included.

The easy way: ``add_organism.py``
---------------------------------

The script ``add_organism.py`` combines the previous two steps (indexing
the fasta files, and adding a gene file), and makes sure the gene BED
file is in the correct place with the correct name. This is the easiest
way to add a new genome/organism for use with GimmeMotifs.


.. _adding_subtools:

Adding motif prediction tools
-----------------------------

Please note that these steps are only necessary when you have installed
any of these tools after you have installed GimmeMotifs.

Weeder
~~~~~~

After installing Weeder the following section needs to be added to the
GimmeMotifs configuration file:

::

    [Weeder]
    bin = /usr/share/Weeder/weederTFBS.out
    dir = /usr/share/Weeder/ 

All other Weeder binaries should be present in the same directory as
``weederTFBS.out``. The directory specified by ``dir`` should contain
the FreqFiles directory included with Weeder. In addition ``Weeder``
should be added to the line in the ``params`` section of the
configuration file. For instance

::

    tools = MDmodule,MEME,MotifSampler,trawler,Improbizer,BioProspector

needs to be changed to:

::

    tools = MDmodule,MEME,MotifSampler,trawler,Improbizer,BioProspector,Weeder

.. _MotifSampler:

MotifSampler configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~

If you want to use MotifSampler there is one more step that you’ll have
to take *after* installation of GimmeMotifs. For every organism, you’ll
need a MotifSampler background. Note that human (hg19), mouse (mm9) and
*Xenopus* (xenTro2) background models are included, so for these
organisms MotifSampler will work out of the box. For other organisms the
necessary background files can be created with ``CreateBackgroundModel``
(which is included with GimmeMotifs or can be downloaded from the same
site as MotifSampler). The background model file needs to be saved in
the directory ``/usr/share/gimmemotifs/MotifSampler`` and it should be
named ``<organism_index_name>.bg``. So, for instance, if I downloaded
the human epd background
(``epd_homo_sapiens_499_chromgenes_non_split_3.bg``), this file should
be saved as ``/usr/share/gimmemotifs/MotifSampler/hg19.bg``.

Other configuration options
---------------------------

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
    background = genomic_matched,random
    use_strand = False
    tools = MDmodule,Weeder,MotifSampler
    analysis = medium
    pvalue = 0.001
    width = 200
    fraction = 0.2
    genome = hg19
    lwidth = 500
    cluster_threshold = 0.95
    available_tools = Weeder,MDmodule,MotifSampler,gadem,meme,trawler
    abs_max = 1000
    enrichment = 1.5
    max_time = None
    scan_cutoff = 0.9

This section specifies all the default GimmeMotifs parameters. Most of
these can also be specified at the command-line when running
GimmeMotifs, in which case they will override the parameters specified
here.
