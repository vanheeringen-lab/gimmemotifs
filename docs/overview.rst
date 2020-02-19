.. _`overview`:

Overview
========

Running on a cluster
--------------------

The current implementation of the cache that GimmeMotifs uses does not play nice with concurrent access, for instance on a cluster. The result is that the cache will get corrupted and that GimmeMotifs will fail. Until this is fixed, there is a workaround. In your job submission script, use something like the following:

::

    NEW_CACHE=$TMPDIR/cache
    mkdir -p $NEW_CACHE
    if [ -z $XDG_CACHE_HOME ]; then
        XDG_CACHE_HOME=$HOME/.cache
    fi
    cp -r $XDG_CACHE_HOME/gimmemotifs $NEW_CACHE/
    export XDG_CACHE_HOME=$NEW_CACHE
    echo 'Using $XDG_CACHE_HOME for cache'

This will use a local directory to store the cache. 

Motif databases
---------------

By default GimmeMotifs uses a non-redundant, clustered database of known vertebrate motifs: **gimme.vertebrate.v5.0**. These motifs come from CIS-BP (http://cisbp.ccbr.utoronto.ca/) and other sources. Large-scale benchmarks using ChIP-seq peaks show that this database shows good performance and should be a good default choice.

In addition, many other motif databases come included with GimmeMotifs:

* `CIS-BP` - All motifs from the `CIS-BP database`_ (version 1.02).
* `ENCODE` - `ENCODE`_ motifs from `Kheradpour & Kellis (2013)`_.
* `factorbook` - All non-redundant motifs from `Factorbook`_.
* `HOCOMOCOv11_HUMAN` - All human motifs from HOCOMOCO_ version 11.
* `HOCOMOCOv11_MOUSE` - All mouse motifs from HOCOMOCO_ version 11.
* `HOMER` - All motifs from HOMER_ (downloaded Oct. 2018).
* `IMAGE` - The motif database from `Madsen et al. (2018)`_.
* `JASPAR2018` - All CORE motifs from `JASPAR 2018`_.
* `JASPAR2018_vertebrates` - CORE vertebrates motifs from `JASPAR 2018`_.
* `JASPAR2018_plants` - CORE plants motifs from `JASPAR 2018`_.
* `JASPAR2018_insects` - CORE insects motifs from `JASPAR 2018`_.
* `JASPAR2018_fungi` -CORE fungi motifs from `JASPAR 2018`_.
* `JASPAR2018_nematodes` - CORE nematodes motifs from `JASPAR 2018`_.
* `JASPAR2018_urochordata` - CORE urochordata motifs from `JASPAR 2018`_.
* `JASPAR2020` - All CORE motifs from `JASPAR 2020`_.
* `JASPAR2020_vertebrates` - CORE vertebrates motifs from `JASPAR 2020`_.
* `JASPAR2020_plants` - CORE plants motifs from `JASPAR 2020`_.
* `JASPAR2020_insects` - CORE insects motifs from `JASPAR 2020`_.
* `JASPAR2020_fungi` -CORE fungi motifs from `JASPAR 2020`_.
* `JASPAR2020_nematodes` - CORE nematodes motifs from `JASPAR 2020`_.
* `JASPAR2020_urochordata` - CORE urochordata motifs from `JASPAR 2020`_.
* `SwissRegulon` - The `SwissRegulon`_ motifs.

You can specify any of these motif databases by name in any GimmeMotifs tool. For instance: 

::

    $ gimme scan TAp73alpha.fa -p JASPAR2018_vertebrates -g hg38

or 

::

    $ gimme motifs TAp73alpha.fa TAp73alpha.motifs -p HOMER -g hg38 --known

.. _`Kheradpour & Kellis (2013)`: https://dx.doi.org/10.1093/nar/gkt1249 
.. _`Madsen et al. (2018)`: https://dx.doi.org/10.1101/gr.227231.117
.. _`Factorbook`: http://www.factorbook.org/human/chipseq/tf/
.. _`ENCODE`: http://compbio.mit.edu/encode-motifs/
.. _`CIS-BP database`: http://cisbp.ccbr.utoronto.ca/
.. _`JASPAR 2018`: http://jaspar.genereg.net
.. _`JASPAR 2020`: http://jaspar.genereg.net
.. _HOMER: http://homer.ucsd.edu/homer/motif/
.. _HOCOMOCO: http://hocomoco11.autosome.ru/
.. _`SwissRegulon`: http://swissregulon.unibas.ch/sr/

