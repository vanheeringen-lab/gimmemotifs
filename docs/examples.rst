.. _`simple_examples`:

Simple examples
===============

Install a genome
----------------

Any genome on UCSC, Ensembl or NCBI can be installed automatically using genomepy_. The ``genomepy`` command tool comes installed with gimmemotifs. For instance, to download the hg38 genome from UCSC:

::

    $ genomepy install hg38 UCSC --annotation


.. _genomepy: https://github.com/simonvh/genomepy

Predict de novo motifs
----------------------

For a quick analysis, to see if it works:

::
    
    $ gimme motifs TAp73alpha.bed -g hg19 -n gimme.p73 -a small -t Homer,MDmodule,BioProspector

For a full analysis:

::
    
    $ gimme motifs TAp73alpha.bed -g hg19 -n gimme.p73 


Open the file ``gimme.p73/gimme.denovo.html`` in your web browser to see the results.

Compare motifs between data sets
--------------------------------

::

    $ gimme maelstrom hg19.blood.most_variable.1k.txt hg19 maelstrom.out/

The output scores of ``gimme maelstrom`` represent the combined result of multiple methods. 
The individual results from different methods are ranked from high-scoring motif to low-scoring motif and converted
to z-scores using the inverse normal transformation. The z-scores from individual methods are then combined using
Stouffer's method. The score that is shown is the aggregated z-score. A higher z-score means that presence of 
the motif or a higher motif score is associated with higher signal in a specific sample.

Create sequence logos
---------------------

::

    $ gimme logo -i GM.5.0.p53.0004

This will generate ``GM.5.0.p53.0004.png``.

.. image:: images/GM.5.0.p53.0004.png

Use the ``-p`` argument for a different pwm file. 
The following command will generate a logo for every motif in ``custom.pfm``.

::  

   $ gimme logo -p custom.pfm

 
