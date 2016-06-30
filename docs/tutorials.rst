Tutorials
=========

Find de novo motifs
-------------------

gimme motifs

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


Installing genomes
------------------

To use most of the functionality of GimmeMotifs you will need to install a genome. 
Is your genome of interest on UCSC? Then you're in luck. If not, don't despair. 
It's still pretty easy, just a few more steps.

Installing from UCSC: ::

    $ gimme genome /data/genomes/ hg38 
    
This will do several things. First, the `FASTADIR` argument, `/data/genomes` in the example above,
determines where the genome FASTA files will be stored. Be aware that this is the genome `root`
directory. A subdirectory with the genome name will be created here.
The second argument specifies the UCSC genome build. 
In this case, `hg38` is the latest version of the human genome on UCSC.
All the genomes that the UCSC Genome Browser supports sohuld be installable in this way. 
Just a few more examples:

Install the Drosophila genome, in a subdir of my home directory: ::

    $ gimme genome ~/genomes sacCer3
    
Install the zebrafish genome, in the current directory: ::

    $ gimme genome . danRer7
    

Installing a non-UCSC genome: 

* Download the FASTA file
* Create a directory with one sequence per file
* gimme index

Create sequence logos
---------------------

gimme logo

Motif enrichment statistics
---------------------------

gimme roc






