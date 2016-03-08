Installation
============

GimmeMotifs runs on Linux. Definitely not on Windows, sorry. Mac OS X
should work in theory, but as I don’t have the means to test this, I’m
not completely sure.

Using pip
---------

Installation from PyPI with `pip` is the most straightforward option:

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
[documentation](https://virtualenv.readthedocs.org/en/latest/userguide.html#usage).

Install numpy:

::

    $ pip install numpy


Now you can install GimmeMotifs using pip. Latest stable release:

::

    $ pip install gimmemotifs


Installation packages
---------------------

Installation packages for Ubuntu and Fedora are no longe supported.

Installation from source
------------------------

Prerequisites
~~~~~~~~~~~~~

Before you can install GimmeMotifs you’ll need:

-  some Python modules and other packages

-  motif prediction tools

Required packages (Python)
~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Python 2.5, 2.6 or 2.7 (not Python 3) http://www.python.org

-  | Scipy http://www.scipy.org/
   | SciPy is the fundamental package needed for scientific computing
   with Python.

-  | matplotlib (0.98 or higher) http://matplotlib.sourceforge.net/
   | A python 2D plotting library. All figures and plots produced by
   GimmeMotifs are made using matplotlib.

-  | kid http://www.kid-templating.org/
   | A simple template language for XML based vocabularies; used to
   produce the HTML reports.

Other required packages
~~~~~~~~~~~~~~~~~~~~~~~

-  ghostscript

Additional motif prediction programs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A lot of motif prediction tools are compiled and/or installed with
GimmeMotifs. The following tools have to be installed seperately:

-  Weeder http://159.149.109.9/modtools/

Please consult the respective manuals regarding installation of these
tools. It’s always possible to install these programs after installation
of GimmeMotifs and update the configuration files to include the new
tools (see section :ref:`adding_subtools`). However, during
installation, GimmeMotifs will try to find any installed tools and add
them automatically, so that’s the easiest option.

Building from source
~~~~~~~~~~~~~~~~~~~~

You can download the lastest stable version of GimmeMotifs at:

| https://github.com/simonvh/gimmemotifs/releases

Start by unpacking the source archive

::

    tar xvzf gimmemotifs-0.8.9.1.tar.gz
    cd gimmemotifs-0.8.9.1

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


