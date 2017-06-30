.. _`api`:

API documentation
=================

.. toctree::
    :maxdepth: 2

    api

Examples
========

Working with motifs
-------------------

The Motif class stores motif information. 
There are several ways to create a Motif instance.

.. code-block:: python

    from gimmemotifs.motif import Motif,read_motifs
    
    # Read from file
    with open("example.pwm") as f:
        motifs = read_motifs(f)
    
    for motif in motifs:
        print(motif)

::

    AP1_nTGAGTCAy
    CTCF_CCAsyAGrkGGCr

.. code-block:: python

    # Create from scratch
    m = Motif([[0,1,0,0],[0,0,1,0]])
    m.id = "CpG"
    print(m)

::

    CpG_CG

.. code-block:: python

    # Or from a consensus sequence
    from gimmemotifs.motif import motif_from_consensus
    ap1 = motif_from_consensus("TGASTCA")
    print(ap1.to_pwm())

::
    
    >TGASTCA
    0.0001  0.0001  0.0001  0.9998
    0.0001  0.0001  0.9998  0.0001
    0.9998  0.0001  0.0001  0.0001
    0.0001  0.4999  0.4999  0.0001
    0.0001  0.0001  0.0001  0.9998
    0.0001  0.9998  0.0001  0.0001
    0.9998  0.0001  0.0001  0.0001

Read motifs from files in other formats.

.. code-block:: python

    with open("motifs.txt") as f:
        motifs = read_motifs(f, fmt="jaspar")

You can convert a motif to several formats.

.. code-block:: python

    with open("example.pwm") as f:
        motifs = read_motifs(f)

    # pwm
    print(motifs[0].to_pwm())

:: 

    >AP1
    0.4908  0.1862  0.2475  0.0755
    0.0125  0.0102  0.0179  0.9594
    0.0191  0.0151  0.9236  0.0422
    0.9457  0.0349  0.0037  0.0158
    0.0355  0.2714  0.6704  0.0228
    0.0121  0.0023  0.0052  0.9804
    0.0271  0.9665  0.0042  0.0022
    0.9935  0.0018  0.0021  0.0027
    0.0367  0.2994  0.1227  0.5412
    
.. code-block:: python

    # pfm
    print(motifs[0].to_pfm())

::

    >AP1
    490.836673106   186.173418152   247.513020751   75.4768879912
    12.547339755300001      10.155349184    17.9452120263   959.3520990339999
    19.1226148238   15.059923645000001      923.609905966   42.2075555647
    945.664342835   34.901043985399994      3.65113750964   15.783475669799998
    35.5141768941   271.353190693   670.372169269   22.7604631439
    12.0584747671   2.26583975448   5.24634274557   980.429342733
    27.0894652078   966.542517084   4.185842951     2.18217475676
    993.460370648   1.82289509211   2.05495304279   2.6617812170899997
    36.6674279259   299.433215935   122.739692952   541.1596631870001

.. code-block:: python

    # consensus sequence
    print(motifs[0].to_consensus())

::

    nTGAGTCAy

Some other useful tidbits.
    
.. code-block:: python

    m = motif_from_consensus("NTGASTCAN")
    print(len(m))

::

    9


.. code-block:: python

    # Trim by information content
    m.trim(0.5)
    print(m.to_consensus(), len(m))

::

    TGAsTCA 7

.. code-block:: python
    
    # Slices
    print(m[:3].to_consensus())
    
::

    TGA

.. code-block:: python
    
    # Shuffle
    random_motif = motif_from_consensus("NTGASTGAN").randomize()
    print(random_motif)

::

    random_snCTAGTAn

To convert a motif to an image, use ``to_img()``.
Supported formats are png, ps and pdf.


.. code-block:: python

    m = motif_from_consensus("NTGASTCAN")
    m.to_img("ap1.png", fmt="png")

.. image:: images/ap1.png

Motif scanning
--------------

For very simple scanning, you can just use a Motif instance. 
Let's say we have a FASTA file called ``test.fa`` that looks like this:

::

    >seq1
    AAAAAAAAAAAAAAAAAAAAAA
    >seq2
    CGCGCGTGAGTCACGCGCGCGCG
    >seq3
    TGASTCAAAAAAAAAATGASTCA

Now we can use this file for scanning.

.. code-block:: python
    
    from gimmemotifs.motif import motif_from_consensus
    from gimmemotifs.fasta import Fasta
    
    f = Fasta("test.fa")
    m = motif_from_consensus("TGAsTCA")

    m.pwm_scan(f)

::

    {'seq1': [], 'seq2': [6, 6], 'seq3': [0, 16, 0, 16]}

This return a dictionary with the sequence names as keys. 
The value is a list with positions where the motif matches. 
Here, as the AP1 motif is a palindrome, you see matches on both forward and reverse strand. 
This is more clear when we use ``pwm_scan_all()`` that returns position, score and strand for every match.

.. code-block:: python

    m.pwm_scan_all(f)

::

    {'seq1': [],
     'seq2': [(6, 9.02922042678255, 1), (6, 9.02922042678255, -1)],
     'seq3': [(0, 8.331251500673487, 1),
     (16, 8.331251500673487, 1),
     (0, 8.331251500673487, -1),
     (16, 8.331251500673487, -1)]}

The number of matches to return is set to 50 by default, you can control this by setting the ``nreport`` argument. 
Use ``scan_rc=False`` to only scan the forward orientation.

.. code-block:: python

    m.pwm_scan_all(f, nreport=1, scan_rc=False)
    
::

    {'seq1': [],
     'seq2': [(6, 9.02922042678255, 1)],
     'seq3': [(0, 8.331251500673487, 1)]}

While this functionality works, it is not very efficient. 
To scan many motifs in potentially many sequences, use the functionality in the ``scanner`` module.
If you only want the best match per sequence, is a utility function called ``scan_to_best_match``, otherwise, use the ``Scanner`` class.

.. code-block:: python

    
    from gimmemotifs.motif import motif_from_consensus
    from gimmemotifs.scanner import scan_to_best_match
    
    m1 = motif_from_consensus("TGAsTCA")
    m1.id = "AP1"
    m2 = motif_from_consensus("CGCG")
    m2.id = "CG"
    motifs = [m1, m2]

    print("motif\tpos\tscore")
    result = scan_to_best_match("test.fa", motifs)
    for motif, matches in result.items():
        for match in matches:
            print("{}\t{}\t{}".format(motif, match[1], match[0]))
    
::

    motif       pos     score
    CG  0       -18.26379789133924
    CG  0       5.554366880674296
    CG  0       -7.743307225501047
    AP1 0       -20.052563923836903
    AP1 6       9.029486018303187
    AP1 0       8.331550321011443

The matches are in the same order as the sequences in the original file.

While this function can be very useful, a ``Scanner`` instance is much more flexible. 
You can scan different input formats (BED, FASTA, regions), and control the thresholds and output.

As an example we will use the file ``Gm12878.CTCF.top500.w200.fa`` that contains 500 top CTCF peaks.
We will get the CTCF motif and scan this file in a number of different ways.

.. code-block:: python

    from gimmemotifs.motif import default_motifs
    from gimmemotifs.scanner import Scanner
    from gimmemotifs.fasta import Fasta
    import numpy as np
    
    # Input file
    fname = "examples/Gm12878.CTCF.top500.w200.fa"
    
    # Select the CTCF motif from the default motif database
    motifs = [m for m in default_motifs() if "CTCF" in m.factors]
    
    # Initialize the scanner
    s = Scanner()
    s.set_motifs(motifs)

Now let's get the best score for the CTCF motif for each sequence.

.. code-block:: python

    scores = [r[0] for r in s.best_score("examples/Gm12878.CTCF.top500.w200.fa")]
    print("{}\t{:.2f}\t{:.2f}\t{:.2f}".format(
        len(scores), 
        np.mean(scores), 
        np.min(scores), 
        np.max(scores)
        ))
    
::

    500	10.61	1.21	14.16

In many cases you'll want to set a threshold. 
In this example we'll use a 1% FPR threshold, based on scanning randomly selected sequences from the ghg38 genome.
The first time you run this, it will take a while. 
However, the tresholds will be cached. 
This means that for the same combination of motifs and genome, the previously generated threshold will be used.

.. code-block:: python

    # Set a 1% FPR threshold based on random hg38 sequence
    s.set_threshold(fpr=0.01, genome="hg38")
    
    # get the number of sequences with at least one match
    counts = [n[0] for n in s.count("examples/Gm12878.CTCF.top500.w200.fa", nreport=1)]
    print(counts[:10])

::

    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

.. code-block:: python

    # or the grand total of number of sequences with 1 match
    print(s.total_count("examples/Gm12878.CTCF.top500.w200.fa", nreport=1))

::

    [408]

.. code-block:: python

    # Scanner.scan() just gives all information
    seqs = Fasta("examples/Gm12878.CTCF.top500.w200.fa")[:10]
    for i,result in enumerate(s.scan(seqs)):
        seqname = seqs.ids[i]
        for m,matches in enumerate(result):
            motif = motifs[m]
            for score, pos, strand in matches:
                print(seqname, motif, score, pos, strand)

::

    chr11:190037-190237 C2H2_ZF_Average_200_CCAsyAGrkGGCr 13.4959558370929 143 -1
    chr11:190037-190237 C2H2_ZF_Average_200_CCAsyAGrkGGCr 10.821440417077262 22 -1
    chr11:190037-190237 C2H2_ZF_Average_200_CCAsyAGrkGGCr 10.658439190070851 82 -1
    chr14:106873577-106873777 C2H2_ZF_Average_200_CCAsyAGrkGGCr 14.16061638444734 120 -1
    chr14:106873577-106873777 C2H2_ZF_Average_200_CCAsyAGrkGGCr 13.72460285196088 83 -1
    chr14:106873577-106873777 C2H2_ZF_Average_200_CCAsyAGrkGGCr 11.450778540447134 27 -1
    chr14:106873577-106873777 C2H2_ZF_Average_200_CCAsyAGrkGGCr 10.037330832055455 7 -1
    chr14:106873577-106873777 C2H2_ZF_Average_200_CCAsyAGrkGGCr 8.998038360035828 159 -1
    chr14:106873577-106873777 C2H2_ZF_Average_200_CCAsyAGrkGGCr 8.668660161058972 101 -1
    chr14:106765204-106765404 C2H2_ZF_Average_200_CCAsyAGrkGGCr 14.16061638444734 145 -1
    chr14:106765204-106765404 C2H2_ZF_Average_200_CCAsyAGrkGGCr 13.848270770440264 185 -1
    chr14:106765204-106765404 C2H2_ZF_Average_200_CCAsyAGrkGGCr 13.668171128367552 165 -1
    chr14:106765204-106765404 C2H2_ZF_Average_200_CCAsyAGrkGGCr 12.785329839873164 27 -1
    chr14:106765204-106765404 C2H2_ZF_Average_200_CCAsyAGrkGGCr 11.886792072933595 126 -1
    chr14:106765204-106765404 C2H2_ZF_Average_200_CCAsyAGrkGGCr 11.25063146496227 67 -1
    chr15:22461178-22461378 C2H2_ZF_Average_200_CCAsyAGrkGGCr 14.16061638444734 28 -1
    chr15:22461178-22461378 C2H2_ZF_Average_200_CCAsyAGrkGGCr 14.16061638444734 185 -1
    chr15:22461178-22461378 C2H2_ZF_Average_200_CCAsyAGrkGGCr 13.261096435278661 67 -1
    chr15:22461178-22461378 C2H2_ZF_Average_200_CCAsyAGrkGGCr 11.450778540447134 147 -1
    chr15:22461178-22461378 C2H2_ZF_Average_200_CCAsyAGrkGGCr 11.022594547749485 126 -1
    chr15:22461178-22461378 C2H2_ZF_Average_200_CCAsyAGrkGGCr 10.194691222675097 7 -1
    chr14:107119996-107120196 C2H2_ZF_Average_200_CCAsyAGrkGGCr 11.886792072933595 37 -1
    chr14:107119996-107120196 C2H2_ZF_Average_200_CCAsyAGrkGGCr 11.886792072933595 95 -1
    chr14:107119996-107120196 C2H2_ZF_Average_200_CCAsyAGrkGGCr 11.886792072933595 153 -1
    chr14:107119996-107120196 C2H2_ZF_Average_200_CCAsyAGrkGGCr 9.972530270193543 75 -1
    chr14:107119996-107120196 C2H2_ZF_Average_200_CCAsyAGrkGGCr 9.949273408029276 17 -1
    chr14:107119996-107120196 C2H2_ZF_Average_200_CCAsyAGrkGGCr 9.949273408029276 133 -1
    chr14:107238917-107239117 C2H2_ZF_Average_200_CCAsyAGrkGGCr 14.16061638444734 92 -1
    chr14:107238917-107239117 C2H2_ZF_Average_200_CCAsyAGrkGGCr 11.25063146496227 34 -1
    chr14:107238917-107239117 C2H2_ZF_Average_200_CCAsyAGrkGGCr 9.246743494108388 15 -1
    chr6:53036754-53036954 C2H2_ZF_Average_200_CCAsyAGrkGGCr 8.764279993851783 62 1
    chr14:107147705-107147905 C2H2_ZF_Average_200_CCAsyAGrkGGCr 13.697109967765122 33 -1
    chr14:107147705-107147905 C2H2_ZF_Average_200_CCAsyAGrkGGCr 13.204664711685334 149 -1
    chr14:107147705-107147905 C2H2_ZF_Average_200_CCAsyAGrkGGCr 11.222131525579154 92 -1
    chr14:107147705-107147905 C2H2_ZF_Average_200_CCAsyAGrkGGCr 11.222131525579154 130 -1
    chr14:50328834-50329034 C2H2_ZF_Average_200_CCAsyAGrkGGCr 11.148765667117496 133 1
    chr1:114889205-114889405 C2H2_ZF_Average_200_CCAsyAGrkGGCr 9.752478102244137 124 1

.. _`api documentation`:

Finding de novo motifs
----------------------

Let's take the ``Gm12878.CTCF.top500.w200.fa`` file as example again. 
For a basic example we'll just use two motif finders, as they're quick to run.

.. code-block:: python

    from gimmemotifs.denovo import gimme_motifs

    peaks = "Gm12878.CTCF.top500.w200.fa"
    outdir = "CTCF.gimme"
    params = {
        "tools": "Homer,BioProspector",
        }

    motifs = gimme_motifs(peaks, outdir, params=params)

::

    2017-06-30 07:37:00,079 - INFO - starting full motif analysis
    2017-06-30 07:37:00,082 - INFO - preparing input (FASTA)
    2017-06-30 07:37:32,949 - INFO - starting motif prediction (medium)
    2017-06-30 07:37:32,949 - INFO - tools: BioProspector, Homer
    2017-06-30 07:37:40,540 - INFO - BioProspector_width_5 finished, found 5 motifs
    2017-06-30 07:37:41,308 - INFO - BioProspector_width_7 finished, found 5 motifs
    2017-06-30 07:37:41,609 - INFO - BioProspector_width_6 finished, found 5 motifs
    2017-06-30 07:37:42,003 - INFO - BioProspector_width_8 finished, found 5 motifs
    2017-06-30 07:37:44,054 - INFO - Homer_width_5 finished, found 5 motifs
    2017-06-30 07:37:45,201 - INFO - Homer_width_6 finished, found 5 motifs
    2017-06-30 07:37:48,246 - INFO - Homer_width_7 finished, found 5 motifs
    2017-06-30 07:37:50,503 - INFO - Homer_width_8 finished, found 5 motifs
    2017-06-30 07:37:54,649 - INFO - BioProspector_width_9 finished, found 5 motifs
    2017-06-30 07:37:56,169 - INFO - BioProspector_width_10 finished, found 5 motifs
    2017-06-30 07:37:56,656 - INFO - Homer_width_9 finished, found 5 motifs
    2017-06-30 07:37:59,313 - INFO - Homer_width_10 finished, found 5 motifs
    2017-06-30 07:37:59,314 - INFO - all jobs submitted
    2017-06-30 07:39:21,298 - INFO - predicted 60 motifs
    2017-06-30 07:39:21,326 - INFO - 53 motifs are significant
    2017-06-30 07:39:21,410 - INFO - clustering significant motifs.
    2017-06-30 07:39:47,031 - INFO - creating reports
    2017-06-30 07:40:41,024 - INFO - finished
    2017-06-30 07:40:41,024 - INFO - output dir: CTCF.gimme
    2017-06-30 07:40:41,024 - INFO - report: CTCF.gimme/motif_report.html

This will basically run the same pipeline as the ``gimme motifs`` command.
All output files will be stored in ``outdir`` and ``gimme_motifs`` returns a list of Motif instances.
If you only need the motifs but not the graphical report, you can decide to skip it by setting ``create_report`` to ``False``.
Additionally, you can choose to skip clustering (``cluster=False``) or to skip calculation of significance (``filter_significant=False``). 
For instance, the following command will only predict motifs and cluster them.

.. code-block:: python

    gimme_motifs(peaks, outdir,
        params=params, filter_significant=False, create_report=False)

All parameters for motif finding are set by the ``params`` argument. 


Motif statistics
----------------

Maelstrom
---------


Auto-generated
==============

Prediction of de novo motifs
----------------------------

.. automodule:: gimmemotifs.denovo
   :members:

Motif scanning
--------------

.. automodule:: gimmemotifs.scanner
   :members:

Motif activity prediction
-------------------------

.. automodule:: gimmemotifs.moap
   :members:

Motif statistics
----------------

.. automodule:: gimmemotifs.stats
   :members: 
