.. _`api`:

API Examples
============

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





.. _`api documentation`:

API Documentation
=================







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
