.. _`overview`:

Overview
========

Motif databases
---------------

By default GimmeMotifs uses a non-redundant, clustered database of known vertebrate motifs: **gimme.vertebrate.v5.0**. These motifs come from CIS-BP (http://cisbp.ccbr.utoronto.ca/) and other sources. 

Several other motif databases come included with GimmeMotifs:

* `CIS-BP` - All motifs from the `CIS-BP database`_ (version 1.02).
* `ENCODE` - 
* `factorbook` - 
* `HOCOMOCOv11_HUMAN` - All human motifs from HOCOMOCO_ version 11.
* `HOCOMOCOv11_MOUSE` - All mouse motifs from HOCOMOCO_ version 11.
* `HOMER` - All motifs from HOMER_ (downloaded Oct. 2018).
* `IMAGE` - 
* `JASPAR2018` - All CORE motifs from `JASPAR 2018`_.
* `JASPAR2018_vertebrates` - CORE vertebrates motifs from `JASPAR 2018`_.
* `JASPAR2018_plants` - CORE plants motifs from `JASPAR 2018`_.
* `JASPAR2018_insects` - CORE insects motifs from `JASPAR 2018`_.
* `JASPAR2018_fungi` -CORE fungi motifs from `JASPAR 2018`_.
* `JASPAR2018_nematodes` - CORE nematodes motifs from `JASPAR 2018`_.
* `JASPAR2018_urochordata` - CORE urochordata motifs from `JASPAR 2018`_.
* `SwissRegulon` - 

You can specify any of these motif databases by name in any GimmeMotifs tool. For instance: 

::

    $ gimme scan TAp73alpha.fa -p JASPAR2018_vertebrates

or 

::

    $ gimme roc TAp73alpha.fa bg.fa -p HOMER -r roc.report


.. _`CIS-BP database`: http://cisbp.ccbr.utoronto.ca/
.. _`JASPAR 2018`: http://jaspar.genereg.net
.. _HOMER: http://homer.ucsd.edu/homer/motif/
.. _HOCOMOCO: http://hocomoco11.autosome.ru/


