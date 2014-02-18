# GimmeMotifs in virtualenv

## Ubuntu

To install GimmeMotifs in a virtualenv, several Python packages need to be built from source. This requires the following packages:

To build Python packages:

* python-dev
* build-essential

To build scipy:

* gfortran
* libblas-dev
* liblapack-dev
* libatlas-base-dev
* cython

To build matplotlib:

* libpng12-dev
* libfreetype6-dev

To build GimmeMotifs:

* libgsl0-dev

Install the necessary packages:

    sudo apt-get install python-pip python-dev build-essential libatlas-base-dev gfortran liblapack-dev libatlas-base-dev cython libpng12-dev libfreetype6-dev libgsl0-dev

Now you can install GimmeMotifs using pip. Latest stable release:

    pip install https://github.com/simonvh/gimmemotifs/tarball/0.8.2

Or the (unstable) master branch with the newest bells, whistles and bugs:

    pip install https://github.com/simonvh/gimmemotifs/tarball/master

