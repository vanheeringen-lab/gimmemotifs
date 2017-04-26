Tutorials
=========

While GimmeMotifs was originally developed to predict *de novo* motifs in ChIP-seq peaks, it is now a full-fledged suite of TF motif analysis tools. 
You can still identify new motifs, but also scan for known motifs, find differential motifs in multiple sets of sequences, create sequence logos, calculate all kinds of enrichment statistics, and more!

Find de novo motifs
-------------------

gimme motifs

Scan for known motifs
---------------------

To scan for known motifs, you will need a set op input sequences and a file with motifs. By default, `gimme scan` uses the motif database that comes included, which is based on clustered, non-redundant motifs from CIS-BP. 
For input sequences you can use either a BED file or a FASTA file. 
You will also need to define the genome. This is used to retrieve sequences, if you have specified a BED file, but also to determine the correct motif-specific threshold for scanning. The default is the genome that is defined in the configuration file.

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






