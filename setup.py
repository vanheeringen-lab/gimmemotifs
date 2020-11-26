from setuptools import setup, find_packages
from setuptools import Extension

import os
from io import open
from compile_externals import compile_all
import versioneer

CONFIG_NAME = "gimmemotifs.cfg"
DESCRIPTION = "GimmeMotifs is a motif prediction pipeline."

with open("README.md", encoding="utf-8") as f:
    long_description = f.read().strip("\n")

# are we in the conda build environment?
conda_build = os.environ.get("CONDA_BUILD")

module1 = Extension("gimmemotifs.c_metrics", sources=["gimmemotifs/c_metrics.c"])

MOTIF_BINS = {
    "MDmodule": ["src/MDmodule/MDmodule"],
    "BioProspector": ["src/BioProspector/BioProspector"],
    "Posmo": ["src/posmo/posmo", "src/posmo/clusterwd"],
    "AMD": ["src/AMD/AMD.bin"],
}


cmdclass = versioneer.get_cmdclass()
my_build_py = cmdclass["build_py"]


class build_tools(my_build_py):
    user_options = []

    def initialize_options(self):
        my_build_py.initialize_options(self)
        self.build_base = None
        self.build_lib = None

    def finalize_options(self):
        my_build_py.finalize_options(self)
        self.set_undefined_options("build", ("build_base", "build_base"))
        self.set_undefined_options("build", ("build_lib", "build_lib"))

    def run(self):
        if not self.dry_run:
            src_dir = os.path.join(self.build_base, "src")
            target_dir = os.path.join(self.build_lib, "gimmemotifs/included_tools")
            # package_data={'gimmemotifs.data':['data/cfg/*']},

            self.copy_tree("src/", src_dir)
            # mkpath is a distutils helper to create directories
            self.mkpath(target_dir)

            compile_all(src_dir=src_dir)

            for exes in MOTIF_BINS.values():
                for exe in exes:
                    if os.path.exists(exe):
                        self.copy_file(exe, target_dir)
                    exe = os.path.join(self.build_base, exe)
                    if os.path.exists(exe):
                        self.copy_file(exe, target_dir)

            self.copy_tree(
                os.path.join(src_dir, "ChIPMunk"), os.path.join(target_dir, "ChIPMunk")
            )
            self.copy_tree(
                os.path.join(src_dir, "HMS"), os.path.join(target_dir, "HMS")
            )
            self.copy_file(
                os.path.join(src_dir, "MotifSampler/MotifSampler_x86_64"),
                os.path.join(target_dir, "MotifSampler"),
            )
            self.copy_file(
                os.path.join(src_dir, "MotifSampler/CreateBackgroundModel_x86_64"),
                os.path.join(target_dir, "CreateBackgroundModel"),
            )
            self.copy_file(
                os.path.join(src_dir, "Improbizer/ameme_x86_64"),
                os.path.join(target_dir, "ameme"),
            )

        my_build_py.run(self)


cmdclass["build_py"] = build_tools

setup(
    name="gimmemotifs",
    version=versioneer.get_version(),
    long_description=long_description,
    long_description_content_type="text/markdown",
    description=DESCRIPTION,
    author="Simon van Heeringen",
    author_email="simon.vanheeringen@gmail.com",
    url="https://github.com/simonvh/gimmemotifs/",
    download_url="https://github.com/simonvh/gimmemotifs/tarball/"
    + versioneer.get_version(),
    license="MIT",
    packages=find_packages(),
    scripts=["scripts/gimme", "scripts/combine_peaks", "scripts/coverage_table"],
    include_package_data=True,
    ext_modules=[module1],
    cmdclass=cmdclass,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.5",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    install_requires=[
        "biofluff",
        "setuptools >= 0.7",
        "numpy",
        "scipy >= 0.9.0",
        "matplotlib >= 2",
        "jinja2",
        "pandas >= 1.1",
        "pyyaml >= 3.10",
        "pybedtools",
        "statsmodels",
        "scikit-learn",
        "seaborn",
        "pysam",
        "xgboost >= 0.71",
        "xdg",
        "diskcache",
        "xxhash",
        "configparser",
        "genomepy >= 0.8.3",
        "tqdm",
        "pillow",
        "logomaker",
        "qnorm",
    ],
)
