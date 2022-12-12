# Copyright (c) 2009-2019 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
"""Interface module for all motif programs."""
import logging
from shutil import which

from .amd import Amd
from .bioprospector import BioProspector
from .chipmunk import ChIPMunk
from .dinamo import Dinamo
from .dreme import Dreme
from .gadem import Gadem
from .hms import Hms
from .homer import Homer
from .improbizer import Improbizer
from .jaspar import Jaspar
from .mdmodule import MDmodule
from .meme import Meme
from .memew import MemeW
from .motifsampler import MotifSampler
from .posmo import Posmo
from .prosampler import ProSampler
from .rpmcmc import Rpmcmc
from .trawler import Trawler
from .weeder import Weeder
from .xxmotif import XXmotif
from .yamda import Yamda

logger = logging.getLogger("gimme.tools")


def get_tool(name):
    """
    Returns an instance of a specific tool.

    Parameters
    ----------
    name : str
        Name of the tool (case-insensitive).

    Returns
    -------
    tool : MotifProgram instance
    """
    tool = name.lower()
    if tool not in __tools__:
        raise ValueError(f"Tool {name} not found!")

    t = __tools__[tool]()

    if not t.is_installed():
        logger.warning(f"Tool {tool} not installed!")

    if not t.is_configured():
        logger.warning(f"Tool {tool} not configured!")

    return t


def locate_tool(name, verbose=True):
    """
    Returns the binary of a tool.

    Parameters
    ----------
    name : str
        Name of the tool (case-insensitive).

    Returns
    -------
    tool_bin : str
        Binary of tool.
    """
    m = get_tool(name)
    tool_bin = which(m.cmd)
    if tool_bin:
        if verbose:
            logger.info(f"Found {m.name} in {tool_bin}")
        return tool_bin
    else:
        if verbose:
            logger.info(f"Couldn't find {m.name}")


__tools__ = {
    "amd": Amd,
    "bioprospector": BioProspector,
    "chipmunk": ChIPMunk,
    "dinamo": Dinamo,
    "dreme": Dreme,
    "gadem": Gadem,
    "homer": Homer,
    "hms": Hms,
    "improbizer": Improbizer,
    "jaspar": Jaspar,
    "mdmodule": MDmodule,
    "meme": Meme,
    "memew": MemeW,
    "motifsampler": MotifSampler,
    "posmo": Posmo,
    "prosampler": ProSampler,
    "rpmcmc": Rpmcmc,
    "trawler": Trawler,
    "weeder": Weeder,
    "xxmotif": XXmotif,
    "yamda": Yamda,
}
