#!/usr/bin/env python
import logging

from yaml import load

from gimmemotifs.tools import get_tool

logger = logging.getLogger("gimme.prediction")


def prediction(args):
    tool = args.tool
    infile = args.infile
    outfile = args.outfile
    paramfile = args.paramfile

    t = get_tool(tool)

    params = {}
    if paramfile:
        params = load(open(paramfile))

    (motifs, stdout, stderr) = t.run(infile, params)

    logger.error(stderr)
    logger.info(stdout)

    f = open(outfile, "w")
    for m in motifs:
        f.write(f"{m.to_pfm()}\n")
    f.close()
