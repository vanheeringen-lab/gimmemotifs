
.. _`command-line`:

Command-line reference
======================

In addition to ``gimme motifs`` the GimmeMotifs package contains
several other tools that can perform the various substeps of
GimmeMotifs, as well as other useful tools. Run them to see the options.

List of tools
-------------

* :ref:`gimme motifs<gimme_motifs>`
* :ref:`gimme maelstrom<gimme_maelstrom>`
* :ref:`gimme scan<gimme_scan>`
* :ref:`gimme roc<gimme_roc>`
* :ref:`gimme match<gimme_match>`
* :ref:`gimme cluster<gimme_cluster>`
* :ref:`gimme index<gimme_index>`
* :ref:`gimme background<gimme_background>`
* :ref:`gimme threshold<gimme_threshold>`
* :ref:`gimme location<gimme_location>`
* :ref:`gimme diff<gimme_diff>`
* :ref:`gimme logo<gimme_logo>`


Input formats
-------------

Most tools in this section take a file in PWM format as input. This is
actually a file with Position Specific Scoring Matrices (PSSMs)
containing *frequencies*. It looks like this:

::

    >motif1
    0.3611  0.0769  0.4003  0.1664
    0.2716  0.0283  0.5667  0.1381
    0.6358  0.0016  0.3344  0.0330
    0.0016  0.9859  0.0016  0.0157
    0.8085  0.0063  0.0502  0.1397
    >motif2
    0.2276  0.0157  0.0330  0.7284
    0.0031  0.0016  0.9984  0.0016
    0.0377  0.3799  0.0016  0.5856
    0.0816  0.7096  0.0173  0.1962
    0.1350  0.4035  0.0675  0.3987

The frequencies are seperated by tabs, and in the order A,C,G,T.


.. _`gimme_motifs`:

Command: gimme motifs
---------------------

Quick example
~~~~~~~~~~~~~

You can try GimmeMotifs with a small example dataset included in the
examples directory, included with GimmeMotifs. This example does not
require any additional configuration if GimmeMotifs is installed
correctly.

Change to a directory where you have write permissions and run the
following command (substitute the filename with the location of the file
on your system):

::

    gimme motifs /usr/share/gimmemotifs/examples/TAp73alpha.fa -n p73

The ``-n`` or ``--name`` option defines the name of the output directory
that is created. All output files are stored in this directory.

Depending on your computer you may have to wait some minutes for your
results. Once GimmeMotifs is finished you can open
`p73/p73\_motif\_report.html <p73/p73_motif_report.html>`__ in your
browser.

Example: gimme motifs
~~~~~~~~~~~~~~~~~~~~~

This example is the same as above, except it will start from a BED file.
This example does require you to have hg19 present and indexed. Change
to a directory where you have write permissions and run the following
command (substitute the filename with the location of the file on your
system):

::

    gimme motifs /usr/share/gimmemotifs/examples/TAp73alpha.bed -n example

The ``-n`` or ``--name`` option defines the name of the output directory
that is created. All output files are stored in this directory.

Depending on your computer you may have to wait some minutes for your
results. Once GimmeMotifs is finished you can open
`example/example\_motif\_report.html <example/example_motif_report.html>`__
in your browser.

Best practices and tips
~~~~~~~~~~~~~~~~~~~~~~~

GimmeMotifs is multi-threaded
+++++++++++++++++++++++++++++

GimmeMotifs runs multi-threaded and uses all the CPU’s in the system.
This means that all the programs will be run in parallel as much as
possible. Of course some programs are still single-threaded, and will
not benefit from this. Because GimmeMotifs uses all the available CPU’s
it does not make much sense to start multiple GimmeMotifs jobs at the
same time.

Running time
++++++++++++

The running time of GimmeMotifs greatly depends on which tools you use
for prediction and how large the dataset is. Some of the tools might
take a very long time and two of them, GADEM is not added to
the default tools because of this reason. You can always use them for an
analysis (by specifying the ``-t`` command-line option), but it is
recommended to only do this for a small dataset (say, less than 5000
peaks). Weeder in combination with the ``xl`` analysis can also take a
very long time, so be prepared. In general a ``small`` analysis will be
the quickest, and a ``xl`` analysis will be the slowest.

While GimmeMotifs is developed specifically for ChIP-seq datasets, most
motif prediction tools are not. In practice this means that it does not
make much sense to predict motifs on a large amount of sequences, as
this will usually not result in higher quality motifs. Therefore
GimmeMotifs uses an absolute limit for the prediction set. By default
20% of the sequences are used as input for motif prediction, but with an
absolute maximum. This is controlled by the ``abs_max`` parameter in the
configuration file, which is set to 1000 by default. In general, if you
have a large amount of peaks, you can also consider to run GimmeMotifs
on the top sequences of your input, for instance the 5000 highest peaks.

There are two options that you can use to control the running time of
GimmeMotifs. First, you can set an absolute time limit with the
``max_time`` option. This option (in hours) determines the maximum time
used for motif prediction. If some programs take longer, the running
jobs will be terminated, and the program will continue with all the
motifs that have been predicted so far. The other option is kind of an
emergency button: when you think that GimmeMotifs has been running long
enough, you can press Ctrl+C **once, and only once!**. This will signal
GimmeMotifs to terminate the running jobs and continue with the
analysis. Please note that this works almost always, but still, there is
a small chance that program might be in a function where the Ctrl-C
option screws up, and GimmeMotifs will not be able to handle the result
gracefully.

Intermediate results
++++++++++++++++++++

GimmeMotifs produces a lot of intermediate results, such as all
predicted motifs, fasta-files used for validation and so on. These are
deleted by default (as they can get quite large), but if you are
interested in them, you can specify the ``-k`` option.

Running on FASTA files
++++++++++++++++++++++

It is possible to run GimmeMotifs on a FASTA file as input instead
of a BED file. This is detected automatically if youir inputfile is
correctly formatted according to FASTA specifications. In this case it
is not possible to generate a genomic matched background, so only the
random Markov background will be used. Please note that for best
results, all the sequences should be of the same length. This is not
necessary for motif prediction, but the statistics and positional
preference plots will be wrong if sequences have different lengths. Also
see the next section.

Small input sets
++++++++++++++++

Keep in mind that GimmeMotifs is developed for larger datasets, where
you have the luxury to use a large fraction of your input for
validation. So, at least several hundred sequences would be optimal. If
you want to run GimmeMotifs on a small input dataset, it might be
worthwile to increase the fraction used for prediction (with the ``-f``
argument.

Detailed options for gimme motifs
+++++++++++++++++++++++++++++++++

-  INPUTFILE

   This is the only mandatory option. The inputfile needs to be in BED
   or FASTA format. BED-fomatted files need to contain at least three
   tab-seperated columns describing chromosome name, start and end. The
   fourth column is optional, if specified it will be used by MDmodule
   to sort the features before motif prediction. GimmeMotifs will take
   the center of these features, and subsequently extend those to the
   width specified by the ``width`` argument (see below).

-  ``-n`` or ``–name``

   The name of your analysis. All outputfiles will be stored in a
   directory named as given by this parameter. By default this will be
   gimmemotifs\_dd\_mm\_yyyy, where d,m and y are the current day, month
   and year respectively.

-  ``-a`` or ``–analysis``

   The size of motifs to look for: small (5-8), medium (5-12), large
   (6-15) or xl (6-20). The larger the motifs, the longer GimmeMotifs
   will run. The ’xl’ can take a very long time!

-  ``-g`` or ``–genome``

   Name of the genome (index) to use. For instance, for the example in
   section :ref:`indexing` this would be ``hg19``.

-  ``-s`` or ``–singlestrand``

   Only use the + strand for prediction (off by default).

-  ``-f`` or ``–fraction``

   This parameter controls the fraction of the sequences used for
   prediction. This 0.2 by default, so in this case a randomly chosen
   20% of the sequences will be used for prediction. The remaining
   sequences will be used for validation (enrichment, ROC curves etc.).
   If you have a large set of sequences (ie. most ChIP-seq peak sets),
   this is fine. However, if your set is smaller, it might be worthwile
   to increase this prediction fraction.

-  ``-w`` or ``–width``

   This is the width of the sequences used for motif prediction. Smaller
   sequences will result in a faster analysis, but you are of course
   limited by the accuracy of your data. For the tested ChIP-seq data
   sets 200 performs fine.

-  ``-e`` or ``–enrichment``

   All motifs should have an absolute enrichment of at least this
   parameter compared to background to be called significant.

-  ``-p`` or ``–pvalue``

   All motifs should have a pvalue of at most this parameter
   (hypergeometric enrichment compared to background) to be called
   significant.

-  ``-b`` or ``–background``

   Type of background to use. By default ``random`` (1st order Markov
   model, similar dinucleotide frequencies as your sequences) and
   ``gc`` (randomly chosen from the genome with a similar
   GC% as your input sequences) are used.

-  ``-l`` or ``–localization_width``

   Width used in the positional preference plots.

-  ``-t`` or ``–tools``

   A comma-seperated list of all the motif prediction tools to use. By
   default all installed tools that are specified in the GimmeMotifs
   configuration file are used.

-  ``–max_time``

   Time limit for motif prediction in hours. Use this to control the
   maximum number of hours that GimmeMotifs uses for motif prediction.
   After this time, all jobs that are still running will be terminated,
   and GimmeMotifs will continue with the motifs that are predicted so
   far.

.. _`gimme_maelstrom`:

Command: gimme maelstrom
------------------------

This command can be used to identify differential motifs between two or more data sets. See the :ref:`maelstrom tutorial<maelstrom_tutorial>` for more details.

**Positional arguments:**

:: 

    INPUTFILE             file with regions and clusters
    GENOME                genome
    DIR                   output directory

**Optional arguments:**

::

    -h, --help            show this help message and exit
    -p PWMFILE, --pwmfile PWMFILE
                          PWM file with motifs (default:
                          gimme.vertebrate.v3.1.pwm)
    -m NAMES, --methods NAMES
                          Run with specific methods

The output scores of `gimme maelstrom` represents the combined result of multiple methods. 
The individual results from different methods are ranked from high-scoring motif to low-scoring motif
and then aggregated using the rank aggregation method from `Kolde, 2012<https://www.ncbi.nlm.nih.gov/pubmed/22247279>`_. 
The score that is shown is the -log10(p-value), where the p-value (from the rank aggregation) is corrected for multiple testing.
This procedure is then repeated with the ranking reversed. These are shown as negative values.

.. _`gimme_scan`:

Command: gimme scan
-------------------

Scan a set of sequences with a set of motifs, and get the resulting
matches in GFF, BED or table format. 
If the FASTA header includes a chromosome location in ``chrom:start-end`` format, the BED output will return the genomic location of the motif match. 
The GFF file will always have the motif location relative to the input sequence.

A basic command would look like this:

::

    $ gimme scan peaks.bed -g hg38 -b > motifs.bed

The threshold that is used for scanning can be specified in a number of ways.
The default threshold is set to a motif-specific 1% FPR by scanning random genomic sequences.
You can change the FPR with the ``-f`` option and/or the set of sequences that is used to determine the FPR with the ``-B`` option.

For instance, this command would scan with thresholds based on 5% FPR with random genomic mouse sequences. 

:: 

    $ gimme scan input.fa -g mm10 -f 0.05 -b > gimme.scan.bed


And this command would base a 0.1% FPR on the input file ``hg38.promoters.fa``:

:: 

    $ gimme scan input.fa -f 0.001 -B hg38.promoters.fa -b > gimme.scan.bed


Alternatively, you can specify the theshold as a single score.
This score is relative and is based on the maximum and minimum possible score for each motif. 
For example, a score of 0.95 means that the score of a motif should be at least 95% of the (maximum score - minimum score).
This should probably not be set much lower than 0.8, and should be generally at least 0.9-0.95 for good specificity. 
Generally, as the optimal threshold might be different for each motif, the use of the FPR-based threshold is preferred.
One reason to use a single score as threshold is when you want a match for each motif, regardless of the score. 
This command would give one match for every motif for every sequence, regardless of the score.

:: 

    $ gimme scan input.bed -g hg38 -c 0 -n 1 -b > matches.bed


Finally, ``gimme scan`` can return the scanning results in table format. 
The ``-t`` will yield a table with number of matches, while the ``-T`` will have the score of the best match.

**Positional arguments:**

:: 

    INPUTFILE             inputfile (FASTA, BED, regions)

**Optional arguments:**

::

    -g GENOME, --genome GENOME
                          genome version
    -p PWMFILE, --pwmfile PWMFILE
                          PWM file with motifs (default:
                          gimme.vertebrate.v3.1.pwm)
    -f , --fpr            FPR for motif scanning (default 0.01)
    -B , --bgfile         background file for threshold
    -c , --cutoff         motif score cutoff or file with cutoffs
    -n N, --nreport N     report the N best matches
    -r, --norc            don't scan reverse complement (- strand)
    -b, --bed             output bed format
    -t, --table           output counts in tabular format
    -T, --score_table     output maximum score in tabular format

.. _`gimme_roc`:

Command: gimme roc
------------------

Given a sample (positives, peaks) and a background file (random
sequences, random promoters or similar), ``gimme roc`` calculates several statistics
and/or creates a ROC plot for motifs in an input PWM file. 
By default, all motifs will be used in the ROC plot, you can select one or more specific motifs with the ``-i`` option. 

The basic command is as follows:

:: 

    $ gimme roc input.fa bg.fa > statistics.txt

This will use the default motif database, and writes the statistics to the file ``statistics.txt``.

Most likely you'll want a graphical report. 
Add the ``-r`` argument to supply an output directory name. 
Once ``gimme roc`` finished, you'll find a file called ``gimme.roc.report.html`` in this directory.
Open it in your browser to get a graphical summary of the results.

Instead of a FASTA file you can also supply a BED file or regions. 
In this case you'll need a genome file.
A custom ``.pwm`` file can be supplied with the ``-p`` argument.
For instance, the following command scans the input BED files with ``custom_motifs.pwm``:

:: 

    $ gimme roc input.bed bg.bed -p custom_motifs.pwm -g hg38 > statistics.txt

The statistics include the ROC area under curve (ROC\_AUC), 
the enrichment at 1% FPR and the recall at 10% FDR.

To plot an ROC curve, add the ``-o`` argument. This command will plot the ROC curve for all the motifs that SPI1 can bind.

::

   $ gimme roc input.fa bg.fa -i Ets_Average_110,Ets_M1778_1.01,Ets_Average_100,Ets_Average_93 -o roc.png > statistics.txt


.. _`Clarke & Granek, 2003`: https://doi.org/10.1093/bioinformatics/19.2.212

**Positional arguments:**

:: 
  
    FG_FILE     FASTA, BED or region file
    BG_FILE     FASTA, BED or region file with background sequences

**Optional arguments:**
  
::

    -h, --help  show this help message and exit
    -r OUTDIR   output dir for graphical report
    -p PWMFILE  PWM file with motifs (default: gimme.vertebrate.v3.1.pwm)
    -g GENOME   Genome (when input files are not in FASTA format)
    -o FILE     Name of output file with ROC plot (png, svg, ps, pdf)
    -i IDS      Comma-seperated list of motif ids to plot in ROC (default is all
                ids)


.. _`gimme_match`:

Command: gimme match
--------------------

Taking an input file with motifs, find the best matching file in another
file of motifs (according to the WIC metric). 
If an ouput file is specified, a graphical output with aligned motifs will
be created. However, this is slow for many motifs and can consume a lot of memory 
(`see issue`_).
It works fine for a few motifs at a time.

.. _`see issue`: https://github.com/simonvh/gimmemotifs/issues/5

**Positional arguments:**

::

    PWMFILE     File with input pwms

**Optional arguments:**

::

    -h, --help  show this help message and exit
    -d DBFILE   File with pwms to match against (default:
                gimme.vertebrate.v3.1.pwm)
    -o FILE     Output file with graphical report (png, svg, ps, pdf)

.. _`gimme_cluster`:

Command: gimme cluster
----------------------

Cluster a set of motifs with the WIC metric.

**Positional arguments:**

::

    INPUTFILE     Inputfile (PFM format)
    OUTDIR        Name of output directory

**Optional arguments:**

::

    -h, --help    show this help message and exit
    -s            Don't compare reverse complements of motifs
    -t THRESHOLD  Cluster threshold

.. _`gimme_index`:

Command: gimme index
--------------------

Creates an index to use with GimmeMotifs.
Use this command if your genome is not available on UCSC and you want to use it with GimmeMotifs.
You should have a directory with FASTA files, **one per chromosome**. 
*Note: this will change with a future version of GimmeMotifs.*

**Positional arguments:**

::

    FASTADIR              Directory to place genome
    GENOMEBUILD           UCSC genome name

**Optional arguments:**

::

    -h, --help            show this help message and exit
    -i DIR, --indexdir DIR
                          Index dir (default
                          <prefix>/share/gimmemotifs/genome_index)


.. _`gimme_background`:

Command: gimme background
-------------------------

Generate random sequences according to one of several methods:

- ``random`` - randomly generated sequence with the same dinucleotide distribution as the input sequences according to a 1st order Markov model
- ``genomic`` - sequences randomly chosen from the genome 
- ``gc`` - sequences randomly chosen from the genome with the same GC% as the input sequences
- ``promoter`` - random promoter sequences

The background types ``gc`` and ``random`` need a set of input sequences
in BED or FASTA format. If the input sequences are in BED format, the 
genome version needs to be specified with ``-g``. 

**Positional arguments:**

::

    FILE        outputfile
    TYPE        type of background sequences to generate
                (random,genomic,gc,promoter)

**Optional arguments:**

::

    -h, --help  show this help message and exit
    -i FILE     input sequences (BED or FASTA)
    -f TYPE     output format (BED or FASTA
    -l INT      length of random sequences
    -n NUMBER   number of sequence to generate
    -g GENOME   genome version (not for type 'random')
    -m N        order of the Markov model (only for type 'random', default 1)

.. _`gimme_threshold`:

Command: gimme threshold
------------------------

Create a file with motif-specific thresholds based on a specific background file and a specific FPR. 
The FPR should be specified as a float between 0.0 and 1.0. 
You can use this threshold file with the ``-c`` argument of :ref:`gimme scan<gimme_scan>`.
Note that :ref:`gimme scan<gimme_scan>` by default determines an FPR based on random genomic background sequences.
You can use this command to create the threshold file explicitly, 
or when you want to determine the threshold based on a different type of background.
For instance, this command would create a file with thresholds for the motifs in ``custom.pwm`` with a FPR of 1%, 
based on the sequences in ``promoters.fa``.

:: 

    $ gimme threshold custom.pwm 0.05 promoters.fa > custom.threshold.txt

**Positional arguments:**

::

    PWMFILE     File with pwms
    FAFILE      FASTA file with background sequences
    FPR         Desired fpr


.. _`gimme_location`:

Command: gimme location
-----------------------

Create the positional preference plots for all the motifs in the input
PWM file. This will give best results if all the sequences in the
FASTA-formatted inputfile have the same length. Keep in mind that this
only makes sense if the sequences are centered around a similar feature
(transcription start site, highest point in a peak, etc.). The default
threshold for motif scanning is 0.95, see ``gimme scan`` for more
details.

**Positional arguments:**

::

    PWMFILE     File with pwms
    FAFILE      Fasta formatted file

**Optional arguments:**

::

    -h, --help  show this help message and exit
    -w WIDTH    Set width to W (default: determined from fastafile)
    -i IDS      Comma-seperated list of motif ids to plot (default is all ids)
    -c CUTOFF   Cutoff for motif scanning (default 0.95)



.. _`gimme_diff`:


Command: gimme diff
-------------------

This is a simple command to visualize differential motifs between different data sets.
You are probably better of using :ref:`gimme maelstrom<gimme_maelstrom>`, however, in some cases this visualization might still be informative.
The input consists of a number of FASTA files, separated by a comma. These are compared to a background file. 
The last two arguments are a file with pwms and and output image. 
The `gimme diff` command then produces two heatmaps (enrichment and frequency) of all enriched, differential motifs.
Reported motifs are at least 3 times enriched compared to the background (change with the ``-e`` argument) and have a minimum frequency in at least one of the input data sets of 1% (change with the ``-f`` argument).
You can specify motif threshold with the ``-c`` argument (which can be a file generated with :ref:`gimme threshold<gimme_threshold>`).

For a command like this...

::

    $ gimme diff VEGT_specific.summit.200.fa,XBRA_specific.summit.200.fa,XEOMES_specific.summit.200.fa random.w200.fa gimme_diff_tbox.png -p tbox.pwm -f 0.01 -c threshold.0.01.txt 

...the output will look like this (based on ChIP-seq peaks of T-box factors from `Gentsch et al. 2013`_):

.. image:: images/gimme_diff_tbox.png

The image layout is not always optimal. 
If you want to customize the image, you can either save it as a ``.svg`` file, or use the numbers that are printed to stdout. 
The columns are in the same order as the image, the row order may be different as these are clustered before plotting.

Note that the results might differ quite a lot depending on the threshold that is chosen! 
Compare for instance an FPR of 1% vs an FPR of 5%.

.. _`Gentsch et al. 2013`: https://doi.org/10.1016/j.celrep.2013.08.012


**Positional arguments:**

::

    FAFILES               FASTA-formatted inputfiles OR a BED file with an
                          identifier in the 4th column, for instance a cluster
                          number.
    BGFAFILE              FASTA-formatted background file
    PNGFILE               outputfile (image)

**Optional arguments:**

::

    -h, --help            show this help message and exit
    -p PWMFILE, --pwmfile PWMFILE
                          PWM file with motifs (default:
                          gimme.vertebrate.v3.1.pwm)
    -c , --cutoff         motif score cutoff or file with cutoffs (default 0.9)
    -e MINENR, --enrichment MINENR
                          minimum enrichment in at least one of the datasets
                          compared to background
    -f MINFREQ, --frequency MINFREQ
                          minimum frequency in at least one of the datasets
    -g VERSION, --genome VERSION
                          Genome version. Only necessary in combination with a
                          BED file with clusters as inputfile.

.. _`gimme_logo`:

Command: gimme logo
-------------------

Convert one or more motifs in a PWM file to a sequence logo.
You can optionally supply a PWM file, otherwise ``gimme logo`` uses the default.
With the ``-i`` option, you can choose one or more motifs to convert.

This will convert all the motifs in ``CTCF.pwm`` to a sequence logo:

:: 

    $ gimme logo -p CTCF.pwm


This will create logos for ``Ets_Average_100`` and ``Ets_Average_109`` from the default database.

:: 

    $ gimme logo -i Ets_Average_100,Ets_Average_109
