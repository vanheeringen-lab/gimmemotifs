===================================================
GimmeMotifs for transcription factor motif analysis
===================================================

What is GimmeMotifs?
--------------------

GimmeMotifs is an analysis framework for transcription factor motif analysis written in Python.
It contains command-line scripts to predict *de novo* motifs, scan for known motifs, identify differential motifs, calculate motif enrichment statistics, plot sequence logos and more.
In addition, all this functionality is available from a Python API.

GimmeMotifs is free and open source research software. 
If you find it useful please cite our paper:

-  Bruse N and van Heeringen SJ, **GimmeMotifs: an analysis framework for transcription factor motif analysis**,
   bioRxiv, 2018. https://www.biorxiv.org/content/10.1101/474403v1.full. `doi: 10.1101/474403 <https://doi.org/10.1101/474403>`_. 

-  van Heeringen SJ and Veenstra GJC, **GimmeMotifs: a de novo motif
   prediction pipeline for ChIP-sequencing experiments**,
   Bioinformatics. 2011 Jan 15;27(2):270-1. `doi: 10.1093/bioinformatics/btq636
   <http://dx.doi.org/10.1093/bioinformatics/btq636>`_.

Getting started
---------------

* The easiest way to :ref:`install<Install GimmeMotifs>` GimmeMotifs is using bioconda_ on Linux or Mac. From version 0.13.0 only Python 3 (>= 3.4) is supported. 
* Have a look at these :ref:`simple examples<simple_examples>` to get a taste of what is possible.
* Check out the more detailed :ref:`tutorials<tutorials>`.
* Full command-line reference can be found :ref:`here<command-line>`.
* There's also an :ref:`API documentation<api>`.

Get help
--------

* First, check the :ref:`FAQ` for common issues.
* The preferred way to get support is through the `Github issues`_ page.
* Finally, you can reach me by mail_ or via twitter_.

.. _`Github issues`: https://github.com/simonvh/gimmemotifs/issues/
.. _bioconda: https://bioconda.github.io/
.. _mail: simon.vanheeringen@gmail.com
.. _twitter: https://twitter.com/svheeringen

Full contents
-------------

.. toctree::
    :maxdepth: 2

    installation
    overview
    examples
    tutorials
    reference
    api
    faq
    acknowledgments
