# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.
""" Configuration for GimmeMotifs """
import configparser
import sysconfig
import xdg
import os
import logging

logger = logging.getLogger("gimme.config")

### CONSTANTS ###
GM_VERSION = "0.11.0"
BG_TYPES = ["random", "genomic", "gc", "promoter"]
FA_VALID_BGS = ["random", "promoter", "gc", "user", "genomic"]
BED_VALID_BGS = ["random", "genomic", "gc", "promoter", "user"]
BG_RANK = {"user":1, "promoter":2, "gc":3, "random":4, "genomic":5}
FASTA_EXT = [".fasta", ".fa", ".fsa"]

CACHE_DIR = os.path.join(xdg.XDG_CACHE_HOME, "gimmemotifs")

class MotifConfig(object):
    """Configuration object for the gimmemotifs module."""
    __shared_state = {}
    prefix = sysconfig.get_config_var("prefix")
    config_dir = "share/gimmemotifs/gimmemotifs.cfg"
    configs = [
        '/home/docs/checkouts/readthedocs.org/user_builds/gimmemotifs/envs/latest/share/gimmemotifs/gimmemotifs.cfg',
        'cfg/gimmemotifs.cfg.example', 
        os.path.join('/usr', config_dir),
        os.path.join(prefix, config_dir), 
        'build/cfg/gimmemotifs.cfg',
        os.path.expanduser('~/.gimmemotifs.cfg')
    ]
    config = None
    TOOL_SECTION = "tools"
    
    def __init__(self, use_config=""):
        self.__dict__ = self.__shared_state
        if use_config:
            self.config = configparser.ConfigParser()
            cfg = self.config.read(use_config)
        elif not self.config:
            self.config = configparser.ConfigParser()
            cfg = self.config.read(self.configs)
            if not cfg:
                raise ValueError("Configuration file not found!")
    def bin(self, program):
        try:
            exe = self.config.get(program, "bin")
        except: 
            raise ValueError("No configuration found for %s" % program)
        return exe
    
    def set_default_params(self, params):
        if not self.config.has_section("params"):
            self.config.add_section("params")

        for k,v in params.items():
            self.config.set("params", k, v)
    
    def get_default_params(self):
        d = dict(self.config.items("params"))
        for k in ["use_strand", "use_cache"]:
            d[k] = self.config.getboolean("params", k)
        return d

    def get_seqlogo(self):
        try:
            exe = self.config.get("main", "seqlogo")
            return exe
        except Exception:
            return None

    def dir(self, program):
        if self.config.has_section(program):
            if self.config.has_option(program, "dir"):
                try: 
                    return self.config.get(program, "dir")
                except Exception:
                    return None
            else:
                return os.path.dirname(self.bin(program))
        else:
            raise ValueError("No configuration found for %s" % program)
    
    def set_program(self, program, d):
        if not self.config.has_section(program):
            self.config.add_section(program)

        for par,value in d.items():
            self.config.set(program, par, value)

    def set_template_dir(self, path):
        if not self.config.has_section("main"):
            self.config.add_section("main")
        self.config.set("main", "template_dir", path)

    def get_template_dir(self):
        return self.config.get("main", "template_dir")

    def set_score_dir(self, path):
        if not self.config.has_section("main"):
            self.config.add_section("main")
        self.config.set("main", "score_dir", path)

    def get_score_dir(self):
        return self.config.get("main", "score_dir")

    def set_seqlogo(self, exe):
        if not self.config.has_section("main"):
            self.config.add_section("main")
        self.config.set("main", "seqlogo", exe)

    def set_index_dir(self, path):
        if not self.config.has_section("main"):
            self.config.add_section("main")
        self.config.set("main", "index_dir", path)

    def get_index_dir(self):
        return self.config.get("main", "index_dir")

    def set_motif_dir(self, path):
        if not self.config.has_section("main"):
            self.config.add_section("main")
        self.config.set("main", "motif_databases", path)

    def get_motif_dir(self):
        return self.config.get("main", "motif_databases")

    def set_gene_dir(self, path):
        if not self.config.has_section("main"):
            self.config.add_section("main")
        self.config.set("main", "gene_dir", path)

    def get_gene_dir(self):
        return self.config.get("main", "gene_dir")
    
    def set_bg_dir(self, path):
        if not self.config.has_section("main"):
            self.config.add_section("main")
        self.config.set("main", "bg", path)

    def get_bg_dir(self):
        return self.config.get("main", "bg")

    def set_tools_dir(self, path):
        if not self.config.has_section("main"):
            self.config.add_section("main")
        self.config.set("main", "tools", path)

    def get_tools_dir(self):
        return self.config.get("main", "tools")

    def is_configured(self, program):
        return self.config.has_section(program)
    
    def save(self):
        self.config.write(open(os.path.expanduser('~/.gimmemotifs.cfg'), "w"))

    def write(self, fo):
        self.config.write(fo)

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
    
    if params["max_time"]:
        try:
            max_time = params["max_time"] = float(params["max_time"])
        except Exception:
            logger.debug("Could not parse max_time value, setting to no limit")
            params["max_time"] = None

        if params["max_time"] > 0:
            logger.debug("Time limit for motif prediction: %0.2f hours", max_time)
            params["max_time"] = 3600 * params["max_time"]
            logger.debug("Max_time in seconds %0.0f", max_time)
        else:
            logger.debug("Invalid time limit for motif prediction, setting to no limit")
            max_time = params["max_time"]
    else:
        logger.debug("No time limit for motif prediction")

    return params

#if __name__ == "__main__":
#    m = MotifConfig()
#    print m.is_configured("meme")
