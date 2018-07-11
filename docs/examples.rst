.. _`simple_examples`:

Simple examples
===============

Install a genome
----------------

Any genome on UCSC, Ensembl or NCBI can be installed automatically using `genomepy<http://github.com/simonvh/genomepy>`_. The `genomepy` command tools comes installed with gimmemotifs. For instance, to download the hg38 genome from UCSC:

::

    $ genomepy install hg38 UCSC --annotation

Predict de novo motifs
----------------------

For a quick analysis, to see if it works:

::
    
    $ gimme motifs TAp73alpha.bed -g hg19 -n gimme.p73 -a small -t Homer,MDmodule,BioProspector

For a full analysis:

::
    
    $ gimme motifs TAp73alpha.bed -g hg19 -n gimme.p73 


Open the file ``gimme.p73/motif_report.html`` in your web browser to see the results.

Compare motifs between data sets
--------------------------------

::

    $ gimme maelstrom hg19.blood.most_variable.1k.txt hg19 maelstrom.out/

The output scores of `gimme maelstrom` represents the combined result of multiple methods. 
The individual results from different methods are ranked from high-scoring motif to low-scoring motif
and then aggregated using the rank aggregation method from `Kolde, 2012<https://www.ncbi.nlm.nih.gov/pubmed/22247279>`_. 
The score that is shown is the -log10(p-value), where the p-value (from the rank aggregation) is corrected for multiple testing. 
This procedure is then repeated with the ranking reversed. These are shown as negative values.

Create sequence logos
---------------------

::

    $ gimme logo -i p53_Average_8

This will generate ``p53_Avarage_8.png``.

.. image:: images/p53_Average_8.png

Use the ``-p`` argument for a different pwm file. 
The following command will generate a logo for every motif in ``custom.pwm``.

::  

   $ gimme logo -p custom.pwm

 
