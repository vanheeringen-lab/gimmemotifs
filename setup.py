from setuptools import setup
from setuptools import Extension, Command
from setuptools.command.install import install

from distutils.command.build import build
from distutils.command.install import INSTALL_SCHEMES
from distutils.util import get_platform
from distutils import log as dlog
import distutils.sysconfig
from subprocess import Popen
from platform import machine
from gimmemotifs.tools import *
from gimmemotifs.config import *
from gimmemotifs.shutils import which
from glob import glob
import os
import sys
import shutil
from stat import ST_MODE
import time
import inspect
from io import open
from compile_externals import compile_all

CONFIG_NAME = "gimmemotifs.cfg" 
DESCRIPTION  = "GimmeMotifs is a motif prediction pipeline."

with open('README.md', encoding='utf-8') as f:
    long_description = f.read().strip("\n")

# are we in the conda build environment?
conda_build = os.environ.get("CONDA_BUILD")

module1 = Extension('gimmemotifs.c_metrics', sources = ['gimmemotifs/c_metrics.c'])

MOTIF_BINS = {
    "MEME": "src/meme_4.6.0/src/meme.bin",
    "MEMEW": "src/meme_4.6.0/src/meme.bin",
    "MDmodule": "src/MDmodule/MDmodule",
    "BioProspector": "src/BioProspector/BioProspector",
    "Posmo": "src/posmo/posmo",
    "AMD": "src/AMD/AMD.bin",
}

class custom_build(build):
    def run(self):
        build.run(self)
        src_dir = os.path.join(self.build_base, "src")
        if not self.dry_run:
            target_dir = os.path.join(self.build_lib, 'gimmemotifs/included_tools')
            
            self.copy_tree("src/", src_dir)
            # mkpath is a distutils helper to create directories
            #self.mkpath(target_dir)

            with open(os.path.join(target_dir, 'myfile.js'), 'w') as f:
                f.write("hello")
            
            print("!!", src_dir)
            compile_all(src_dir=src_dir)

            for exe in MOTIF_BINS.values():
                if os.path.exists(exe):
                    self.copy_file(exe, target_dir)
            
            self.copy_tree(os.path.join(src_dir,  "ChIPMunk"), os.path.join(target_dir, "ChIPMunk"))
            self.copy_tree(os.path.join(src_dir,"HMS"), os.path.join(target_dir, "HMS"))
            self.copy_file(os.path.join(src_dir,"MotifSampler/MotifSampler_x86_64"), os.path.join(target_dir, "MotifSampler"))
            self.copy_file(os.path.join(src_dir,"MotifSampler/CreateBackgroundModel_x86_64"), os.path.join(target_dir, "CreateBackgroundModel"))
            self.copy_file(os.path.join(src_dir,"Improbizer/ameme_x86_64"), os.path.join(target_dir, "ameme"))

setup (
        name = 'gimmemotifs',
        version = GM_VERSION,
        long_description = long_description,
        long_description_content_type = 'text/markdown',
        description = DESCRIPTION,
        author = 'Simon van Heeringen',
        author_email = 'simon.vanheeringen@gmail.com',
        url = 'https://github.com/simonvh/gimmemotifs/',
        download_url = 'https://github.com/simonvh/gimmemotifs/tarball/' + GM_VERSION,
        license = 'MIT',
        packages=['gimmemotifs', 'gimmemotifs/commands', 'gimmemotifs/included_tools'],
        ext_modules = [module1],
        cmdclass = {
            "build":custom_build,
            },
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            #'License :: OSI Approved :: MIT License',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: POSIX :: Linux',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3.5',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            ],
        scripts=[
            'scripts/track2fasta.py',
            'scripts/gimme',
            ],
        install_requires = [
            "setuptools >= 0.7",
            "numpy",
            "scipy >= 0.9.0",
            "matplotlib >= 2",
            "jinja2",
            "pyyaml >= 3.10",
            "pybedtools",
            "statsmodels",
            "scikit-learn",
            "sklearn-contrib-lightning",
            "seaborn",
            "pysam",
            "xgboost >= 0.71",
            "xdg",
            "diskcache",
            "xxhash",
            "configparser",
            "six",
            "future",
            "genomepy",
            "tqdm",
            ],
)
