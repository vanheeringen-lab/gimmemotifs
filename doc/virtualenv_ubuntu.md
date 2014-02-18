# GimmeMotifs in virtualenv

## Ubuntu

To install GimmeMotifs in a virtualenv, several Python packages need to be built from source. Several packages are necessary:

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

Now you can install GimmeMotifs using pip:

    pip install https://github.com/simonvh/gimmemotifs/tarball/master 



<!---
sudo apt-get install libamd2.2.0 libblas3gf libc6 libgcc1 libgfortran3 liblapack3gf libumfpack5.4.0 libstdc++6 build-essential gfortran libatlas-base-dev python-all-dev

Python
* python-pip
* python-dev
* build-essential




pip install numpy
pip install scipy
pip install kid
pip install matplotlib

pip install pyyaml
