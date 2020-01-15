# Copyright (c) 2009-2019 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
"""Interface module for all motif programs."""
from __future__ import print_function

import sys

from ..shutils import which
from .mdmodule import MDmodule
from .meme import Meme
from .memew import MemeW
from .dreme import Dreme
from .weeder import Weeder
from .gadem import Gadem
from .motifsampler import MotifSampler
from .trawler import Trawler
from .improbizer import Improbizer
from .bioprospector import BioProspector
from .posmo import Posmo
from .chipmunk import ChIPMunk
from .jaspar import Jaspar
from .amd import Amd
from .hms import Hms
from .homer import Homer
from .xxmotif import XXmotif
from .prosampler import ProSampler
from .yamda import Yamda
from .dinamo import Dinamo
from .rpmcmc import Rpmcmc


MOTIF_CLASSES = [
    "MDmodule",
    "Meme",
    "MemeW",
    "Dreme",
    "Weeder",
    "Gadem",
    "MotifSampler",
    "Trawler",
    "Improbizer",
    "BioProspector",
    "Posmo",
    "ChIPMunk",
    "Jaspar",
    "Amd",
    "Hms",
    "Homer",
    "XXmotif",
    "ProSampler",
    "YAMDA",
    "DiNAMO",
    "RPMCMC",
]


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
        raise ValueError("Tool {0} not found!\n".format(name))

    t = __tools__[tool]()

    if not t.is_installed():
        sys.stderr.write("Tool {0} not installed!\n".format(tool))

    if not t.is_configured():
        sys.stderr.write("Tool {0} not configured!\n".format(tool))

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
            print("Found {} in {}".format(m.name, tool_bin))
        return tool_bin
    else:
        print("Couldn't find {}".format(m.name))


__tools__ = {
    "xxmotif": XXmotif,
    "homer": Homer,
    "bioprospector": BioProspector,
    "hms": Hms,
    "amd": Amd,
    "improbizer": Improbizer,
    "trawler": Trawler,
    "weeder": Weeder,
    "motifsampler": MotifSampler,
    "mdmodule": MDmodule,
    "chipmunk": ChIPMunk,
    "posmo": Posmo,
    "gadem": Gadem,
    "jaspar": Jaspar,
    "meme": Meme,
    "memew": MemeW,
    "dreme": Dreme,
    "prosampler": ProSampler,
    "yamda": Yamda,
    "dinamo": Dinamo,
    "rpmcmc": Rpmcmc,
}
