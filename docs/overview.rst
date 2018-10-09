.. _`overview`:

Overview
========


Motif databases
---------------

By default GimmeMotifs uses a non-redundant, clustered database of known vertebrate motifs: **gimme.vertebrate.v3.1**. These motifs come from CIS-BP(http://cisbp.ccbr.utoronto.ca/), and are based on different sources. 

Several other motif databases come included with GimmeMotifs:

* `JASPAR2018` - All CORE motifs from `JASPAR 2018`_.
* `JASPAR2018_vertebrates` - CORE vertebrates motifs from `JASPAR 2018`_.
* `JASPAR2018_plants` - CORE plants motifs from `JASPAR 2018`_.
* `JASPAR2018_insects` - CORE insects motifs from `JASPAR 2018`_.
* `JASPAR2018_fungi` -CORE fungi motifs from `JASPAR 2018`_.
* `JASPAR2018_nematodes` - CORE nematodes motifs from `JASPAR 2018`_.
* `JASPAR2018_urochordata` - CORE urochordata motifs from `JASPAR 2018`_.
* `HOMER` - All motifs from HOMER_ (downloaded Oct. 2018)
* `HOCOMOCOv11_HUMAN` - All human motifs from HOCOMOCO_ version 11.
* `HOCOMOCOv11_MOUSE` - All mouse motifs from HOCOMOCO_ version 11.

You can specify any of these motif databases by name in any GimmeMotifs tool. For instance: 

::

    $ gimme scan TAp73alpha.fa -p JASPAR2018_vertebrates

.. _`JASPAR 2018`: http://jaspar.genereg.net
.. _HOMER: http://homer.ucsd.edu/homer/motif/
.. _HOCOMOCO: http://hocomoco11.autosome.ru/


