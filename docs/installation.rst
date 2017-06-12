Installation
============

GimmeMotifs runs on Linux. Definitely not on Windows, sorry. 
Mac OSX should work and is included in the build test. 
However, as I don't use it myself, unexpected issues might pop up. 
Let me know, so I can try to fix it.

.. _`Install GimmeMotifs`:

The easiest way to install
--------------------------

The most straightforward way to install GimmeMotifs is by using `conda
<https://docs.continuum.io/anaconda>`_. 
Activate the bioconda_ channel if you haven't done so already.

:: 

    $ conda config --add channels conda-forge
    $ conda config --add channels defaults
    $ conda config --add channels r
    $ conda config --add channels bioconda

Now you can install GimmeMotifs with one command. In the current environment:

::

    $ conda install gimmemotifs

Or create a specific environment:

::

    $ conda create -n gimme gimmemotifs
    # Activate the environment before you use GimmeMotifs
    $ source activate gimme

.. _bioconda: https://bioconda.github.io/

Using pip
---------

Installation from PyPI with `pip` is another straightforward option. 
You will need to install some prerequisites though:

- bedtools http://bedtools.readthedocs.io
- UCSC genePredToBed http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToBed
- (optional) R + RobustRankAggreg

Install with pip as follows:

:: 

    $ sudo pip install gimmemotifs

Or the (unstable) master branch with the newest bells, whistles and bugs:

::

    $ sudo pip install https://github.com/simonvh/gimmemotifs/tarball/master

If you don't have root access, see the option below.

Using pip in a virtualenv
-----------------------------

Ubuntu prerequisites
~~~~~~~~~~~~~~~~~~~~

To install GimmeMotifs in a virtualenv, several Python packages need to be built from source. 

Install the necessary packages to build numpy, scipy, matplotlib and GimmeMotifs:

::

    sudo apt-get install python-pip python-dev build-essential libatlas-base-dev \
    gfortran liblapack-dev libatlas-base-dev cython libpng12-dev libfreetype6-dev \
    libgsl0-dev

Install via pip
~~~~~~~~~~~~~~~

Create a virtualenv and activate it according to the 
`documentation
<https://virtualenv.readthedocs.org/en/latest/userguide.html#usage>`_.

Install numpy:

::

    $ pip install numpy


Now you can install GimmeMotifs using pip. Latest stable release:

::

    $ pip install gimmemotifs


Installation packages
---------------------

Installation packages for Ubuntu and Fedora are no longer supported.

Installation from source
------------------------

These instructions are not up-to-date! Basically, you're on your own!

Make sure to install all required dependencies.

You can download the lastest stable version of GimmeMotifs at:

| https://github.com/simonvh/gimmemotifs/releases

Start by unpacking the source archive

::

    tar xvzf gimmemotifs-0.9.0.2.tar.gz
    cd gimmemotifs-0.9.0.2

You can build GimmeMotifs with the following command:

::

    python setup.py build

Run the tests to check if the basics work correctly:

::

    python run_tests.py

If you encounter no errors, go ahead with installing GimmeMotifs (root
privileges required):

::

    sudo python setup.py install

During installation GimmeMotifs will try to locate the tools you have
installed. If you have recently installed them, running an ``updatedb``
will be necessary. Using this option GimmeMotifs will create a
configuration file, the default is:

::

    /usr/share/gimmemotifs/gimmemotifs.cfg

This is a system-wide configuration that can be used by all users.

It is also possible to run the ``setup.py install`` command with the
``–prefix``, ``–home``, or ``–install-data`` options, to install in
GimmeMotifs in a different location (for instance, in your own home
directory). This should be fine, however, these alternative methods of
installing GimmeMotifs have not been extensively tested. Please note
that in this case the configuration file will be created, but every user
will have to put this configuration file in his/her home directory:
``~/.gimmemotifs.cfg``. The install script will also inform you of this
during install. Please contact me if you run into problems with the
installation. Once the installation is finished, you can try the quick
example (section :ref:`quick-example`), or continue with the
configuration in the next section.


