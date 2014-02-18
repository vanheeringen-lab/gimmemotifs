# GimmeMotifs in virtualenv

## Ubuntu

To install GimmeMotifs in a virtualenv, several Python packages need to be built from source. 

Install the necessary packages to build numpy, scipy, matplotlib and GimmeMotifs:

    sudo apt-get install python-pip python-dev build-essential libatlas-base-dev \
    gfortran liblapack-dev libatlas-base-dev cython libpng12-dev libfreetype6-dev \
    libgsl0-dev

Now you can install GimmeMotifs using pip. Latest stable release:

    pip install https://github.com/simonvh/gimmemotifs/tarball/0.8.2

Or the (unstable) master branch with the newest bells, whistles and bugs:

    pip install https://github.com/simonvh/gimmemotifs/tarball/master

