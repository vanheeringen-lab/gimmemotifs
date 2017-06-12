.. _tutorials:

Tutorials
=========

While GimmeMotifs was originally developed to predict *de novo* motifs in ChIP-seq peaks, it is now a full-fledged suite of TF motif analysis tools. 
You can still identify new motifs, but also scan for known motifs, find differential motifs in multiple sets of sequences, create sequence logos, calculate all kinds of enrichment statistics, and more!

For this tutorial I'll assume you use bioconda. 
Create an environment with all the necessary dependencies.

:: 

    $ conda create -n gimme_tutorial gimmemotifs ucsc-bigbedtobed

And activate it!

:: 
    
    $ source activate gimme_tutorial

To locate the example files mentioned in the tutorial, locate the ``examples/`` directory of your GimmeMotifs installation. When using conda:

::

    $  echo `conda info | grep default | awk '{print $4}'`/share/gimmemotifs/examples
    /home/simon/anaconda3/share/gimmemotifs/examples
    $ ls /home/simon/anaconda3/share/gimmemotifs/examples
    TAp73alpha.bed  TAp73alpha.fa

Find de novo motifs
-------------------

As a simple example, let's predict the CTCF motif based on ChIP-seq data from ENCODE.
Download the peaks:

::    

    $ wget http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/integration_data_jan2011/byDataType/peaks/jan2011/spp/optimal/hub/spp.optimal.wgEncodeBroadHistoneGm12878CtcfStdAlnRep0_VS_wgEncodeBroadHistoneGm12878ControlStdAlnRep0.bb

Convert the bigBed file to a BED file using ``bigBedToBed``:

::

    $ bigBedToBed spp.optimal.wgEncodeBroadHistoneGm12878CtcfStdAlnRep0_VS_wgEncodeBroadHistoneGm12878ControlStdAlnRep0.bb Gm12878.CTCF.bed

Select the top 500 peaks, based on the signalValue column of the narrowPeak_ format, as input:

::

    $ sort -k7gr Gm12878.CTCF.bed | awk '{print $1 "\t" $2 + $10 "\t" $2 + $10}' | head -n 500 > Gm12878.CTCF.top500.bed

In addition, we base the location of the peak on the summit (column 10). 
Note the top 500 peaks are just for the tutorial. 
Normally you would use a much larger sample (or all peaks) as input for ``gimme motifs``.

Now, the ENCODE peak coordinates are based on hg19 so we need to install the hg19 genome.
For a UCSC genome, this is just a matter of running ``gimme genome``.

:: 
    
    $ gimme genome /data/genomes/ hg19

I store my genomes in ``/data/genomes`` but you can use any directory you like. 
Also see the section on :ref:`installing genomes<Installing genomes>` below.
This will take some time. 
The genome sequence will be downloaded and indexed, ready for use with GimmeMotifs.

Having both an index genome and an input file, we can run ``gimme motifs``.

:: 

    $ gimme motifs Gm12878.CTCF.top500.bed -g hg19 -n gimme.CTCF

Once again, this will take some time. 
When ``gimme motifs``  is finished you can view the results in a web browser. 
`gimme.CTCF/motif_report.html`_ should look a lot like this.
This is what an almost perfect motif looks like, with a ROC AUC close to 1.0.

.. _`narrowPeak`: https://genome.ucsc.edu/FAQ/FAQformat.html#format12


Scan for known motifs
---------------------

To scan for known motifs, you will need a set of input sequences and a file with motifs. 
By default, ``gimme scan`` uses the motif database that comes included, which is based on clustered, non-redundant motifs from CIS-BP. 
For input sequences you can use either a BED file, a FASTA file or a file with regions in ``chr:start-end`` format. 
You will also need to specify the genome, see the section `Installing genomes`_ below. 
The genome sequence will be used to retrieve sequences, if you have specified a BED or region file, but also to determine a reasonable motif-specific threshold for scanning. 
The default genome can be specified in the configuration file.

We will use the file ``Gm12878.CTCF.top500.bed`` that was used for `de novo` motif search above for known motifs.
While ``gimme motifs`` automatically extends regions from the center of the input regions, ``gimme scan`` uses the regions as specified in the file. 
This means we will have to change the size of the regions to 200 nucleotides. 
Depending on the type and quality of your input data, you can of course make this smaller or larger.

:: 

    $ bedtools slop -i Gm12878.CTCF.top500.bed -b 100 -g hg19.genome > Gm12878.CTCF.top500.w200.bed

OK, let's scan:

::

    $ gimme scan Gm12878.CTCF.top500.w200.bed -g hg19 > result.scan.gff

The first time you run ``gimme scan`` for a specific combination of motif database, genome, input sequence length and FPR (which is 0.01 by default) it will determine a motif-specific cutoff based on random genome background sequences. 
This will take a while. However, results will be cached for future scanning.

To get a BED file with the genomic location of motif matches add the ``-b`` argument:

::

    $ gimme scan Gm12878.CTCF.top500.w200.bed -g hg19 -b > result.scan.bed

By default, ``gimme scan`` gives at most one match per sequence for each motif, if the score of the match reaches a certain threshold.

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

In this case, the most abundant motif is the CTCF motif. 

The specified false positive rate (FPR), with a default of 0.01, determines the motif-specific threshold that is used for scanning.
This means that the expected rate of occurrence, determined by scanning random genomic sequences, is 1%. 
Based on the FPR, you can assume that any motif with more than 1% matches is enriched. 
However, for a more robust measure of motif significance use ``gimme roc``, which is further explained :ref:`below<roc>`.
This command will give the enrichment, but also the ROC AUC and recall at 10% FDR and other useful statistics.

For many applications, it is useful to have motif occurrences as a table. 

:: 

    $ gimme scan Gm12878.CTCF.top500.w200.bed -g hg19 -t > table.count.txt
 
This will result in a tab-separated table with counts. 
Same defaults as above, at most one match per sequence per motif.
Alternatively, ``gimme scan`` can report the score of best match, regardless of the value of this score.

:: 

    $ gimme scan Gm12878.CTCF.top500.w200.bed -g hg19 -T > table.score.txt
    $ head table.score.txt | cut -f1-10
    # GimmeMotifs version 0.10.1b2
    # Input: Gm12878.CTCF.top500.w200.bed
    # Motifs: /home/simon/anaconda3/share/gimmemotifs/motif_databases/gimme.vertebrate.v3.1.pwm
    # FPR: 0.01 (hg19)
            AP-2_Average_26 AP-2_Average_17 AP-2_Average_27 AP-2_Average_15 AP-2_M5965_1.01 ARID_BRIGHT_Average_1   ARID_BRIGHT_M0104_1.01  ARID_BRIGHT_Average_3   ARID_BRIGHT_M5966_1.01
    chr11:190037-190237     3.315682        5.251773        5.852259        6.986044        -0.032952       -1.058302       -4.384525       1.989879        -13.872373
    chr14:106873577-106873777       3.485541        5.315545        3.867055        1.129976        4.386212        -3.305211       -1.392656       2.726421        -13.660561
    chr14:106765204-106765404       3.936576        5.315545        3.867055        1.434064        -1.284617       -1.058302       -3.578581       1.597828        -8.376869
    chr15:22461178-22461378 3.936576        5.315545        3.867055        1.387997        -1.284617       -3.305211       -7.331101       1.551285        -8.929093
    chr14:107119996-107120196       3.485541        5.468490        3.867055        1.434064        4.708942        -5.675314       -7.331101       1.159831        -15.790964

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

.. _roc:

Motif enrichment statistics
---------------------------

You can use ``gimme roc`` to compare motifs or to identify relevant known motifs for a specific input file.

Let's get some known motifs for one of the example files, ``TAp73alpha.fa``. 
First, we need to define a background.
To get random genomic sequences with a matched GC% content:

:: 

    $ gimme background random.gc.fa gc -g hg19 -n 500 -l 200 -i TAp73alpha.fa

This will create a FASTA file with 500 sequences of 200 nucleotides, that has a GC% distribution similar to ``TAp73alpha.fa``.
Now we can run ``gimme roc``:

:: 

    $ gimme roc TAp73alpha.fa random.gc.fa > TAp73alpha.roc.txt

What's in the file?

:: 

    $ head -n 1 TAp73alpha.roc.txt | tr \\t \\n
    Motif
    ROC AUC
    MNCP
    Enr. at 5% FPR
    Max enr.
    Recall at 10% FDR

The motif ID, followed by five statistics: the ROC area under curve (AUC), Mean Normalized Conditional Probability (MNCP), the enrichment compared to background set at 5% FPR, the maximum enrichment and the recall at 10% FDR.

The ROC AUC is widely used, however, it might not always be the most informative.
In situations where the background set is very large compared to the input set, it might give a more optimistic picture than warranted.

Let's sort on the last statistic:

:: 

    $ sort -k6g TAp73alpha.roc.txt | cut -f1,2,6 | tail
    p53_M5923_1.01      0.666   0.1050
    bZIP_M0304_1.01     0.621   0.1240
    Grainyhead_Average_6        0.698   0.1720
    Unknown_M6235_1.01  0.682   0.2330
    bZIP_Average_149    0.620   0.2400
    Myb_SANT_Average_7  0.592   0.2560
    Runt_Average_9      0.709   0.3700
    p53_M3568_1.01      0.823   0.6270
    p53_Average_10      0.825   0.6490
    p53_Average_8       0.918   0.8860

Not surprisingly, the p53 family motif is the most enriched. 
In addition, we also get RUNX1 and AP1 motifs. 
The Grainyhead motif somewhat resembles the p53 motif, which could explain the enrichment. 
Let's visualize this.
This command will create two sequence logos in PNG format:

:: 

    $ gimme logo -i p53_Average_8,Grainyhead_Average_6

The p53 motif, or p73 motif in this case, ``p53_Average_8.png``:

.. image:: images/p53_Average_8.png

And the Grainyhead motif, ``Grainyhead_Average_6``:

.. image:: images/Grainyhead_Average_6.png

The resemblance is clear. 
This also serves as a warning to never take the results from a computational tool (including mine) at face value...
