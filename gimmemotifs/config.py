# Copyright (c) 2009-2019 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
""" Configuration for GimmeMotifs """
import configparser
import glob
import sys
import xdg
import os
import pathlib
import logging
import pkg_resources
from shutil import which
from time import time
from ._version import get_versions

logger = logging.getLogger("gimme.config")

# CONSTANTS #
BG_TYPES = ["random", "genomic", "gc", "promoter"]
FA_VALID_BGS = ["random", "promoter", "gc", "custom", "genomic"]
BED_VALID_BGS = ["random", "genomic", "gc", "promoter", "custom"]
BG_RANK = {"custom": 1, "promoter": 2, "gc": 3, "random": 4, "genomic": 5}
FASTA_EXT = [".fasta", ".fa", ".fsa", ".fna"]
DIRECT_NAME = "direct"
INDIRECT_NAME = "indirect\nor predicted"

CACHE_DIR = os.path.join(xdg.XDG_CACHE_HOME, "gimmemotifs")
CONFIG_DIR = os.path.join(xdg.XDG_CONFIG_HOME, "gimmemotifs")

MOTIF_CLASSES = [
    "AMD",
    "BioProspector",
    "ChIPMunk",
    "DiNAMO",
    "DREME",
    "GADEM",
    "HMS",
    "Homer",
    "Improbizer",
    "MDmodule",
    "MEME",
    "MEMEW",
    "MotifSampler",
    "Posmo",
    "ProSampler",
    "RPMCMC",
    "Trawler",
    "Weeder",
    "XXmotif",
    "Yamda",
]


def get_build_dir():
    """
    Returns the build directory if installed in editable mode
    using `python setup.py build && pip install -e .`

    Returns None if installed regularly using `pip install .`
    """
    root_dir = os.path.dirname(os.path.dirname(__file__))
    v = sys.version_info
    glob_dir = os.path.join(root_dir, "build", f"lib*{v[0]}*{v[1]}*", "gimmemotifs")
    results = glob.glob(glob_dir)

    if len(results) == 1:
        return results[0]


class MotifConfig(object):
    """Configuration object for the gimmemotifs module."""

    # Borg design pattern: all instances of this class will have the same attributes
    __shared_state = {}

    # Default config that is installed with GimmeMotifs
    default_config = pkg_resources.resource_filename(
        "gimmemotifs", "../data/cfg/gimmemotifs.default.cfg"
    )
    user_config = os.path.join(CONFIG_DIR, "gimmemotifs.cfg")
    config = None

    # If gimme is installed in editable mode,
    # the motif discovery tools are installed in the build/ dir,
    # else they are installed in environment's site-packages/ dir.
    package_dir = get_build_dir()
    if package_dir is None:
        package_dir = os.path.dirname(__file__)

    # minimum time before updating the config again
    __timeout = 10  # seconds
    __last_checked = 0

    def __init__(self, use_config=None):
        self.__dict__ = self.__shared_state
        self.config = configparser.ConfigParser()
        if use_config:
            cfg = self.config.read(use_config)
        else:
            cfg = self.config.read(self.user_config)

        if not cfg:
            self.create_default_config()
            self.__last_checked = time()
        elif time() > self.__last_checked + self.__timeout:
            self._upgrade_config()
            self.__last_checked = time()

    def _upgrade_config(self, config_fname=None):
        changed = False
        dflt = configparser.ConfigParser()
        dflt.read(self.default_config)

        # check if old tools are still available
        available_tools = self.config["params"]["available_tools"].split(",")
        for m in available_tools:
            cmd = self.config.get(m, "bin")
            if not os.path.isfile(cmd):
                logger.info("%s no longer found.", m)
                available_tools.remove(m)
                self.set_program(
                    m, {"bin": dflt.get(m, "bin"), "dir": dflt.get(m, "dir")}
                )
                changed = True

        # check if new tools are available
        missing_tools = [m for m in MOTIF_CLASSES if m not in available_tools]
        for m in missing_tools:
            cmd = self.bin(m, config=dflt, missing_ok=True)
            msg = "Using included version of %s."
            if cmd is None:
                cmd = which(dflt.get(m, "bin"))
                msg = "Using system version of %s."
            if cmd:
                logger.info(msg, m)
                available_tools.append(m)
                self.set_program(m, {"bin": cmd, "dir": os.path.dirname(cmd)})
                changed = True

        # update older configs
        if "width" in self.config["params"]:
            if "size" not in self.config["params"]:
                self.config.set(
                    option="size",
                    section="params",
                    value=self.config["params"]["width"],
                )
            del self.config["params"]["width"]
            changed = True
        if "lwidth" in self.config["params"]:
            if "lsize" not in self.config["params"]:
                self.config.set(
                    option="lsize",
                    section="params",
                    value=self.config["params"]["lwidth"],
                )
            del self.config["params"]["lwidth"]
            changed = True

        if changed:
            available_tools.sort()
            self.config["params"]["available_tools"] = ",".join(available_tools)
            if config_fname is None or not config_fname:
                config_fname = self.user_config
            with open(config_fname, "w") as f:
                self.write(f)
            logger.info("Configuration file: %s", self.user_config)

    def create_default_config(self):
        logger.info("Creating new config.")

        available_tools = []
        self.config.read(self.default_config)
        self.config.set("main", "config_version", get_versions()["version"])
        for m in MOTIF_CLASSES:
            mbin = self.config.get(m, "bin")
            mdir = self.config.get(m, "dir")
            cmd = which(os.path.join(self.package_dir, mdir, mbin))
            msg = "Using included version of %s."
            if cmd is None:
                cmd = which(mbin)
                msg = "Using system version of %s."
            if cmd:
                logger.info(msg, m)
                self.set_program(m, {"bin": cmd, "dir": os.path.dirname(cmd)})
                available_tools.append(m)
                continue
            logger.warning(
                "%s not found. To include it you will have to install it.", m
            )

        params = self.get_default_params()
        params["available_tools"] = ",".join(available_tools)
        self.set_default_params(params)

        if not os.path.exists(CONFIG_DIR):
            os.makedirs(CONFIG_DIR, exist_ok=True)
        with open(self.user_config, "w") as f:
            self.config.write(f)
        logger.info("Configuration file: %s", self.user_config)

    def bin(self, program, config=None, missing_ok=False):
        if config is None:
            config = self.config

        mbin = config.get(program, "bin")
        if os.path.exists(mbin):
            return os.path.abspath(mbin)

        mdir = config.get(program, "dir")
        cmd = os.path.join(mdir, mbin)
        if os.path.exists(cmd):
            return os.path.abspath(cmd)

        cmd = os.path.join(self.package_dir, mdir, mbin)
        if os.path.exists(cmd):
            return cmd

        if not missing_ok:
            raise ValueError(f"No configuration found for {program}")

    def dir(self, program):
        mdir = self.config.get(program, "dir", fallback="included_tools")
        if os.path.exists(mdir):
            return mdir

        mdir = os.path.join(self.package_dir, mdir)
        if os.path.exists(mdir):
            return mdir

        mdir = os.path.dirname(self.bin(program))
        if os.path.isabs(mdir):
            return mdir

        raise ValueError("No configuration found for %s" % program)

    def set_program(self, program, d):
        if not self.config.has_section(program):
            self.config.add_section(program)

        for par, value in d.items():
            self.config.set(program, par, value)

    def get_data_dir(self, ddir):
        my_dir = self.config.get("main", ddir)
        if not os.path.exists(my_dir):
            my_dir = pkg_resources.resource_filename(
                "gimmemotifs", os.path.join("../data", my_dir)
            )
        return my_dir

    def set_default_params(self, params):
        if not self.config.has_section("params"):
            self.config.add_section("params")

        for k, v in params.items():
            self.config.set("params", k, str(v))

    def get_default_params(self):
        d = dict(self.config.items("params"))
        for k in ["use_strand", "use_cache"]:
            d[k] = self.config.getboolean("params", k)
        return d

    def set_template_dir(self, path):
        if not self.config.has_section("main"):
            self.config.add_section("main")
        self.config.set("main", "template_dir", path)

    def get_template_dir(self):
        return self.get_data_dir("template_dir")

    def set_score_dir(self, path):
        if not self.config.has_section("main"):
            self.config.add_section("main")
        self.config.set("main", "score_dir", path)

    def get_score_dir(self):
        return self.get_data_dir("score_dir")

    def set_motif_dir(self, path):
        if not self.config.has_section("main"):
            self.config.add_section("main")
        self.config.set("main", "motif_databases", path)

    def get_motif_dir(self):
        return self.get_data_dir("motif_databases")

    def set_gene_dir(self, path):
        if not self.config.has_section("main"):
            self.config.add_section("main")
        self.config.set("main", "gene_dir", path)

    def get_gene_dir(self):
        return self.get_data_dir("gene_dir")

    def set_bg_dir(self, path):
        if not self.config.has_section("main"):
            self.config.add_section("main")
        self.config.set("main", "bg", path)

    def get_bg_dir(self):
        return self.get_data_dir("bg")

    def set_tools_dir(self, path):
        if not self.config.has_section("main"):
            self.config.add_section("main")
        self.config.set("main", "tools", path)

    def get_tools_dir(self):
        return self.config.get("main", "tools")

    def is_configured(self, program):
        return self.config.has_section(program)

    def write(self, fo):
        self.config.write(fo)

    def list_installed_libraries(self):
        """Return a list of all motif libraries installed in this distribution.

        Each returned string is suitable for use with `gimmemotifs.motif.read_motifs()`.
        """
        libraries_dir = pathlib.Path(self.get_motif_dir()).resolve()
        library_paths = glob.glob(str(libraries_dir / "*.pfm"))
        return sorted([pathlib.Path(p).name for p in library_paths])


def parse_denovo_params(user_params=None):
    """Return default GimmeMotifs parameters.

    Defaults will be replaced with parameters defined in user_params.

    Parameters
    ----------
    user_params : dict, optional
        User-defined parameters.

    Returns
    -------
    params : dict
    """
    config = MotifConfig()

    if user_params is None:
        user_params = {}
    params = config.get_default_params()
    params.update(user_params)

    if params.get("torque"):
        logger.debug("Using torque")
    else:
        logger.debug("Using multiprocessing")

    params["background"] = [x.strip() for x in params["background"].split(",")]

    logger.debug("Parameters:")
    for param, value in params.items():
        logger.debug("  %s: %s", param, value)

    # Maximum time?
    try:
        params["max_time"] = float(params.get("max_time", -1))
    except (ValueError, TypeError):
        logger.debug("Could not parse max_time value, setting to no limit")
        params["max_time"] = -1

    max_time = params["max_time"]
    if max_time > 0:
        logger.debug("Time limit for motif prediction: %0.2f hours", max_time)
        params["max_time"] = 3600 * params["max_time"]
        logger.debug("Max_time in seconds %0.0f", max_time)
    else:
        logger.debug("No time limit for motif prediction")

    return params


# if __name__ == "__main__":
#    m = MotifConfig()
#    print m.is_configured("meme")
