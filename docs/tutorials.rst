Tutorials
=========

Scan for known motifs
---------------------

To scan for known motifs, for instance in a set of ChIP-seq peaks, you will need to do three things:

* Define a genomic background file
* Set a threshold, based on the background
* Scan!

In  short:

Define background::

    $ gimme background hg38.random.w200.fa genomic -n 10000 -w 200 -g hg38

Set threshold (1% FDR)::

    $ gimme threshold <pwmfile> background hg38.random.w200.fa 0.01 > hg38.threshold.w200.fdr0.01.txt
    
Scan::

    $ gimme scan <input> -p <pwmfile> -c hg38.threshold.w200.fdr0.01.txt > result.scan.gff

Find differential motifs
------------------------

gimme maelstrom


Compare two sets with de novo motifs
------------------------------------

gimme motifs

combine: gimme cluster

scan
