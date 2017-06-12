.. _tutorials:

Tutorials
=========

While GimmeMotifs was originally developed to predict *de novo* motifs in ChIP-seq peaks, it is now a full-fledged suite of TF motif analysis tools. 
You can still identify new motifs, but also scan for known motifs, find differential motifs in multiple sets of sequences, create sequence logos, calculate all kinds of enrichment statistics, and more!

Throughout this tutorial I'll assume you'll use bioconda.

Find de novo motifs
-------------------

As a simple example, let's predict the CTCF motif based on ChIP-seq data from ENCODE.
Download the peaks:

::    

    $ wget http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/integration_data_jan2011/byDataType/peaks/jan2011/spp/optimal/hub/spp.optimal.wgEncodeBroadHistoneGm12878CtcfStdAlnRep0_VS_wgEncodeBroadHistoneGm12878ControlStdAlnRep0.bb

Install bigBedToBed if you have not done so already:

:: 

    $ conda install ucsc-bigbedtobed

Convert the bigBed file to a BED file:

::

    $ bigBedToBed spp.optimal.wgEncodeBroadHistoneGm12878CtcfStdAlnRep0_VS_wgEncodeBroadHistoneGm12878ControlStdAlnRep0.bb Gm12878.CTCF.bed

Select the top 500 peaks, based on the signalValue column of the narrowPeak format, as input:

::

    $ sort -k7gr Gm12878.CTCF.bed | awk '{print $1 "\t" $2 + $10 "\t" $2 + $10}' | head -n 500 > Gm12878.CTCF.top500.bed

In addition, we base the location of the peak on the summit (column 10). 
Note the top 500 peaks are just for the tutorial. 
Normally you would use a much larger sample or all peaks as input for ``gimme motifs``.

Now, the ENCODE peaks are based on hg19 so we need to install the hg19 genome.
For a UCSC genome, this is just a matter of running ``gimme genome``.

:: 
    
    $ gimme genome /data/genomes/ hg19

I store my genomes in ``/data/genomes`` but you can use any directory you like. 
Also see the section `Installing genomes`_ below.
This will take some time. 
The genome sequence will be downloaded and indexed, ready for use with GimmeMotifs.

Having both an index genome and an input file, we can run ``gimme motifs``.

:: 

    $ gimme motifs Gm12878.CTCF.top500.bed -g hg19


Scan for known motifs
---------------------

To scan for known motifs, you will need a set op input sequences and a file with motifs. 
By default, ``gimme scan`` uses the motif database that comes included, which is based on clustered, non-redundant motifs from CIS-BP. 
For input sequences you can use either a BED file, a FASTA file or a file with regions in ``chr:start-end`` format. 
You will also need to specify the genome, see the section `Installing genomes`_ below. 
The genome sequence will be used to retrieve sequences, if you have specified a BED or region file, but also to determine the correct motif-specific threshold for scanning. 
The default is the genome that is defined in the configuration file.

We will use the file ``Gm12878.CTCF.top500.bed`` that was used for `de novo` motif search above for known motifs.
Change the size of the regions to 200 nucleotides:

:: 

    $ bedtools slop -i Gm12878.CTCF.top500.bed -b 100 -g hg19.genome > Gm12878.CTCF.top500.w200.bed

Scan: 

::

    $ gimme scan Gm12878.CTCF.top500.w200.bed -g hg19 > result.scan.gff

The first time you run ``gimme scan`` for a specific combination of motif database, genome, input sequence length and FPR (which is 0.01 by default) it will determine a motif-specific cutoff based on random genome background sequences. 
This will take a while. However, results will be cached for future scanning.

To get a BED file with the genomic location of motif matches add the ``-b`` argument:

::

    $ gimme scan Gm12878.CTCF.top500.w200.bed -g hg19 -b > result.scan.bed

For a very simple summary, we can just have a look at the most abundant motifs:

:: 

    $ cut -f4 result.scan.bed | sort | uniq -c | sort -n | tail
         50 E2F_Average_31
         58 C2H2_ZF_Average_123
         58 MBD_Average_1
         58 THAP_finger_M1541_1.01
         59 Unknown_Average_5
         67 Ets_Average_70
         72 C2H2_ZF_M0401_1.01
         72 CxxC_M0548_1.01
        118 E2F_Average_27
        394 C2H2_ZF_Average_200

The most abundant motif is the CTCF motif. 


Find differential motifs
------------------------

gimme maelstrom


Compare two sets with de novo motifs
------------------------------------

gimme motifs

combine: gimme cluster

scan

.. _`Installing genomes`:

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






