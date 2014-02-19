Installation
============

GimmeMotifs runs on Linux. Definitely not on Windows, sorry. Mac OS X
should work in theory, but as I don’t have the means to test this, I’m
not completely sure. I have tried to make installation of GimmeMotifs as
easy as possible. There are installation packages available for Fedora,
Ubuntu and Debian. It’s also possible to install GimmeMotifs from
source, but it depends on quite some external packages. Please make sure
all prerequisites are installed before installing GimmeMotifs from
source.

Installation packages
---------------------

Installation on Ubuntu or Debian
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Download the latest GimmeMotifs .deb package for your distribution from
http://www.ncmls.nl/bioinfo/gimmemotifs/ (available for 32-bit and
64-bit systems). Please note: due to a different Python version, the
Debian and Ubuntu .deb packages are different! Please use the correct
one for your distribution. These packages have been tested on up-to-date
versions of Ubuntu Oneiric Ocelot (11.10) and Debian sqeeze (6.04). Open
a terminal and install GimmeMotifs as follows. Go to the directory where
you downloaded GimmeMotifs, for example:

::

    cd ~/Downloads

Start the install (substitute the package name for the .deb package you
downloaded):

::

    sudo dpkg -i gimmemotifs_0.65_amd64.deb 

Likely, dpkg will complain about some missing dependencies. Install all
dependencies with:

::

    sudo apt-get -f install

Complete the GimmeMotifs installation with:

::

    sudo dpkg -i gimmemotifs_0.65_amd64.deb 

Currently, there is a bug with the versions of the Parallel Python
(python-pp) and Numpy (python-numpy) in the Ubuntu and Debian
repositories. Therefore the package python-pp is not installed as a
dependency. This can be fixed by installing the latest version of
Parallel Python from the Python Package Index:

::

    sudo easy_install pp 

Now you should have a working version of GimmeMotifs! The next steps are
to install additional motif tools (optional, see section
:ref:`adding_subtools`) and to do some configuration (required, see
section :ref:`Configuration` ). You can also directly try the quick
example (:ref:`quick-example`), if you’re impatient (but
don’t forget to perform the additional steps!)

Installation on Fedora
~~~~~~~~~~~~~~~~~~~~~~

Download the latest GimmeMotifs .rpm package from
http://www.ncmls.nl/bioinfo/gimmemotifs/ (available for 32-bit and
64-bit systems). This package has been tested on an up-to-date version
of Fedora 16. Install GimmeMotifs as follows (substitute the package
name for the .rpm package you downloaded):

::

    sudo yum install --nogpgcheck gimmemotifs-0.65-1.x86_64.rpm 

GimmeMotifs doesn’t play nice with SELinux enabled on Fedora, sorry. To
turn it off:

::

    sudo setenforce 0

Now you should have a working version of GimmeMotifs. The next steps are
to install additional motif tools (optional, see section
:ref:`adding_subtools`) and to do some configuration (required, see
:ref:`Configuration` ). You can also directly try the quick
example (see :ref:`quick-example`), if you’re impatient (but
don’t forget to perform the additional steps!)

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

-  | parallel python 1.6.0 http://www.parallelpython.com/
   | A python module which provides mechanism for parallel execution of
   python code. This Python library is used for parallel execution of
   for instance the motif finding tools.

-  | kid http://www.kid-templating.org/
   | A simple template language for XML based vocabularies; used to
   produce the HTML reports.

Other required packages
~~~~~~~~~~~~~~~~~~~~~~~

-  | gsl http://www.gnu.org/software/gsl/
   | The GNU Scientific Library. This library might already be installed
   on your system, but you’ll also need the development headers to
   compile GimmeMotifs!.

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

| You can download the lastest version of GimmeMotifs at:
| http://www.ncmls.eu/bioinfo/gimmemotifs/.
| Start by unpacking the source archive

::

    tar xvzf gimmemotifs-1.00.tar.gz
    cd gimmemotifs-1.00

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
