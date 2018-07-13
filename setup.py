from setuptools import setup, find_packages
from setuptools import Extension, Command
from distutils.command.build import build
from setuptools.command.install import install

import os
import glob
import sys
from io import open
from compile_externals import compile_all

GM_VERSION = "0.12.1"

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

class build_tools(Command):
    def initialize_options(self):
        self.build_base = None
        self.build_lib = None

    def finalize_options(self):
        self.set_undefined_options('build',('build_base', 'build_base'))
        self.set_undefined_options('build',('build_lib', 'build_lib'))

    def run(self): 
        if not self.dry_run:
            src_dir = os.path.join(self.build_base, "src")
            target_dir = os.path.join(self.build_lib, 'gimmemotifs/included_tools')
            
            self.copy_tree("src/", src_dir)
            # mkpath is a distutils helper to create directories
            self.mkpath(target_dir)
            
            compile_all(src_dir=src_dir)

            for exe in MOTIF_BINS.values():
                if os.path.exists(exe):
                    self.copy_file(exe, target_dir)
            
            self.copy_tree(
                    os.path.join(src_dir,  "ChIPMunk"), 
                    os.path.join(target_dir, "ChIPMunk"))
            self.copy_tree(
                    os.path.join(src_dir,"HMS"), 
                    os.path.join(target_dir, "HMS"))
            self.copy_file(
                    os.path.join(src_dir,"MotifSampler/MotifSampler_x86_64"), 
                    os.path.join(target_dir, "MotifSampler"))
            self.copy_file(
                    os.path.join(src_dir,"MotifSampler/CreateBackgroundModel_x86_64"), 
                    os.path.join(target_dir, "CreateBackgroundModel"))
            self.copy_file(
                    os.path.join(src_dir,"Improbizer/ameme_x86_64"), 
                    os.path.join(target_dir, "ameme"))
       
            if os.path.exists("src/weblogo"):
                self.copy_tree("src/weblogo", 
                        os.path.join(target_dir, "weblogo"))

build.sub_commands += [
            ('build_tools', lambda self: True),
            ]

class custom_install(install):
 
    def run(self):
        self.run_command('build_tools')
        self.do_egg_install()

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
        packages=find_packages(),
        scripts=[
            'scripts/track2fasta.py',
            'scripts/gimme',
            ],
        #package_data={'gimmemotifs.data':['data/cfg/*']},
        include_package_data = True,
        ext_modules = [module1],
        cmdclass = {
            "install":custom_install,
            "build_tools":build_tools,
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
