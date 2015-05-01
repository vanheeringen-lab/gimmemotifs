Usage
=====

.. _quick-example:

Quick GimmeMotifs example
-------------------------

You can try GimmeMotifs with a small example dataset included in the
examples directory, included with GimmeMotifs. This example does not
require any additional configuration if GimmeMotifs is installed
correctly.

Change to a directory where you have write permissions and run the
following command (substitute the filename with the location of the file
on your system):

::

    gimme motifs /usr/share/gimmemotifs/examples/TAp73alpha.fa -n p73

The ``-n`` or ``–name`` option defines the name of the output directory
that is created. All output files are stored in this directory.

Depending on your computer you may have to wait some minutes for your
results. Once GimmeMotifs is finished you can open
`p73/p73\_motif\_report.html <p73/p73_motif_report.html>`__ in your
browser.

GimmeMotifs example
-------------------

This example is the same as above, except it will start from a BED file.
This example does require you to have hg19 present and indexed. Change
to a directory where you have write permissions and run the following
command (substitute the filename with the location of the file on your
system):

::

    gimme motifs /usr/share/gimmemotifs/examples/TAp73alpha.bed -n example

The ``-n`` or ``–name`` option defines the name of the output directory
that is created. All output files are stored in this directory.

Depending on your computer you may have to wait some minutes for your
results. Once GimmeMotifs is finished you can open
`example/example\_motif\_report.html <example/example_motif_report.html>`__
in your browser.

Using GimmeMotifs: best practices and tips
------------------------------------------

GimmeMotifs is multi-threaded
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GimmeMotifs runs multi-threaded and uses all the CPU’s in the system.
This means that all the programs will be run in parallel as much as
possible. Of course some programs are still single-threaded, and will
not benefit from this. Because GimmeMotifs uses all the available CPU’s
it does not make much sense to start multiple GimmeMotifs jobs at the
same time.

Running time
~~~~~~~~~~~~

The running time of GimmeMotifs greatly depends on which tools you use
for prediction and how large the dataset is. Some of the tools might
take a very long time and two of them, GADEM and MoAn, are not added to
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
~~~~~~~~~~~~~~~~~~~~

GimmeMotifs produces a lot of intermediate results, such as all
predicted motifs, fasta-files used for validation and so on. These are
deleted by default (as they can get quite large), but if you are
interested in them, you can specify the ``-k`` option.

Running on FASTA files
~~~~~~~~~~~~~~~~~~~~~~

It is also possible to run GimmeMotifs on a FASTA file as input instead
of a BED file. This is detected automatically if you’re inputfile
correctly formatted according to FASTA specifications. In this case it
is not possible to generate a genomic matched background, so only the
random Markov background will be used. Please note that for best
results, all the sequences should be of the same length. This is not
necessary for motif prediction, but the statistics and positional
preference plots will be wrong if sequences have different lengths. Also
see the next section.

Small input sets
~~~~~~~~~~~~~~~~

Keep in mind that GimmeMotifs is developed for larger datasets, where
you have the luxury to use a large fraction of your input for
validation. So, at least several hundred sequences would be optimal. If
you want to run GimmeMotifs on a small input dataset, it might be
worthwile to increase the fraction used for prediction (with the ``-f``
parameter.

Detailed options
----------------

-  INPUTFILE

   This is the only mandatory option. The inputfile needs to be in BED
   or FASTA format. BED-fomatted files need to contain at least three
   tab-seperated columns describing chromosome name, start and end. The
   fourth column is optional, if specified it will be used by MDmodule
   to sort the features before motif prediction. GimmeMotifs will take
   the center of these features, and subsequently extend those to the
   width specified by the ``width`` parameter (see below).

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
