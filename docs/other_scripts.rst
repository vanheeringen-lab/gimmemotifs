Other scripts
=============

In addition to ``gimme motifs`` the GimmeMotifs package contains
several other tools that can perform the various substeps of
GimmeMotifs, as well as other useful tools. Run them to see the options.

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

Descriptions
------------

gimme match
~~~~~~~~~~~

Taking an input file with motifs, find the best matching file in another
file of motifs (according to the WIC metric).

gimme index
~~~~~~~~~~~

Creates an index to use with GimmeMotifs. See :ref:`Configuration` for details.

gimme background
~~~~~~~~~~~~~~~~

Generate random sequences according to one of two methods: random or
matched\_genomic. With the argument ``type`` set to ``random``, and an
input file in FASTA format, this script will generate sequences with the
same dinucleotide distribution as the input sequences according to a 1st
order Markov model trained on the input sequences. The ``-n`` options is
set to 10 by default. The length distribution of the sequences in the
output file will be similar as the inputfile. The Markov model can be
changed with option ``-m``. If the ``type`` is specified as
``matched_genomic`` the inputfile needs to be in BED format, and the
script will select genomic regions with a similar distribution relative
to the transcription start of genes as the input file. Make sure to
select the correct genome. The length of the sequences in the output
file will be set to the median of the features in the input file.

gimme cluster
~~~~~~~~~~~~~

Cluster a set of motifs with the WIC metric.

gimme location
~~~~~~~~~~~~~~

Create the positional preference plots for all the motifs in the input
PWM file. This will give best results if all the sequences in the
FASTA-formatted inputfile have the same length. Keep in mind that this
only makes sense if the sequences are centered around a similar feature
(transcription start site, highest point in a peak, etc.). The default
threshold for motif scanning is 0.95, see ``gimme scan`` for more
details.

gimme roc
~~~~~~~~~

Given a sample (positives, peaks) and a background file (random
sequences, random promoters or similar), calculates several statistics
and/or creates a ROC plot for all the motifs in an input PWM file. All
the motifs will be plotted in the same graph, you can select one or more
specific motifs to plot with the ``-i`` option. The statistics include
ROC area under curve (ROC\_AUC) and Mean Normalized Conditional
Probability (MNCP).

gimme logo
~~~~~~~~~~

Convert the motifs in a PWM file to a logo using weblogo.

gimme scan
~~~~~~~~~~

Scan a set of sequences with a set of motifs, and give the resulting
matches in GFF or BED format. The threshold is based on the maximum and
minimum possible score for each motif. So, 0.95 means that the score of
a motif should be at least 95% of the (maximum score - minimum score).
This should probably not be set much lower than 0.8, and should be
generally at least 0.9 for good specificity. Keep in mind that the
optimal threshold might be different for each motif!

track2fasta.py 
~~~~~~~~~~~~~~~

Convert a set of BED formatted sequences to a FASTA file. The genome
needs to be indexed for GimmeMotifs using ``gimme index``.
