import logging
import os
import sys

from genomepy import Genome

from gimmemotifs.config import BED_VALID_BGS, FA_VALID_BGS
from gimmemotifs.utils import determine_file_type

logger = logging.getLogger()


def check_bed_file(fname):
    """Check if the inputfile is a valid bed-file"""
    if not os.path.exists(fname):
        logger.error(f"Inputfile {fname} does not exist!")
        sys.exit(1)

    for i, line in enumerate(open(fname)):
        if (
            line.startswith("#")
            or line.startswith("track")
            or line.startswith("browser")
        ):
            # comment or BED specific stuff
            pass
        else:
            vals = line.strip().split("\t")
            if len(vals) < 3:
                logger.error(
                    "Expecting tab-seperated values (chromosome<tab>start<tab>end) "
                    f"on line {i + 1} of file {fname}",
                )
                sys.exit(1)
            try:
                int(vals[1]), int(vals[2])
            except ValueError:
                logger.error(
                    f"No valid integer coordinates on line {i + 1} of file {fname}"
                )
                sys.exit(1)
            if len(vals) > 3:
                try:
                    float(vals[3])
                except ValueError:
                    pass


def check_denovo_input(inputfile, params):
    """
    Check if an input file is valid, which means BED, narrowPeak or FASTA
    """
    background = params["background"]

    input_type = determine_file_type(inputfile)

    if input_type == "fasta":
        valid_bg = FA_VALID_BGS
    elif input_type in ["bed", "narrowpeak"]:
        genome = params["genome"]
        valid_bg = BED_VALID_BGS
        if "genomic" in background or "gc" in background:
            Genome(genome)
        # is it a valid bed-file etc.
        check_bed_file(inputfile)  # bed-specific, will also work for narrowPeak
    else:
        logger.error(f"Format of inputfile {inputfile} not recognized.")
        logger.error("Input should be FASTA, BED or narrowPeak.")
        logger.error(
            "See https://genome.ucsc.edu/FAQ/FAQformat.html for specifications."
        )
        sys.exit(1)

    for bg in background:
        if bg not in valid_bg:
            logger.info(f"Input type is {input_type}, ignoring background type '{bg}'")
        background = [bg for bg in background if bg in valid_bg]

    if len(background) == 0:
        logger.error("No valid backgrounds specified!")
        sys.exit(1)

    return input_type, background
