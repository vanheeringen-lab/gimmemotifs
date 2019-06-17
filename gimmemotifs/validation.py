from gimmemotifs.utils import determine_file_type
from gimmemotifs.config import FA_VALID_BGS, BED_VALID_BGS
from genomepy import Genome
import os
import sys
import logging

logger = logging.getLogger()


# import logger


def check_bed_file(fname):
    """ Check if the inputfile is a valid bed-file """
    if not os.path.exists(fname):
        logger.error("Inputfile %s does not exist!", fname)
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
                    "Expecting tab-seperated values (chromosome<tab>start<tab>end) on line %s of file %s",
                    i + 1,
                    fname,
                )
                sys.exit(1)
            try:
                int(vals[1]), int(vals[2])
            except ValueError:
                logger.error(
                    "No valid integer coordinates on line %s of file %s", i + 1, fname
                )
                sys.exit(1)
            if len(vals) > 3:
                try:
                    float(vals[3])
                except ValueError:
                    pass
                    # self.logger.warn("No numerical value in column 4 on line %s of file %s, ignoring..." % (i + 1, file))


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
        sys.stderr.write("Format of inputfile {} not recognized.\n".format(inputfile))
        sys.stderr.write("Input should be FASTA, BED or narrowPeak.\n")
        sys.stderr.write(
            "See https://genome.ucsc.edu/FAQ/FAQformat.html for specifications.\n"
        )
        sys.exit(1)

    for bg in background:
        if bg not in valid_bg:
            logger.info(
                "Input type is %s, ignoring background type '%s'", input_type, bg
            )
        background = [bg for bg in background if bg in valid_bg]

    if len(background) == 0:
        logger.error("No valid backgrounds specified!")
        sys.exit(1)

    return input_type, background
