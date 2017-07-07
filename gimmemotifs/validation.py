import os
import sys
import logging

logger = logging.getLogger()


from gimmemotifs.fasta import Fasta
from gimmemotifs.genome_index import check_genome
from gimmemotifs.config import (MotifConfig, FA_VALID_BGS, BED_VALID_BGS)
# import logger

def check_bed_file(fname):
    """ Check if the inputfile is a valid bed-file """
    if not os.path.exists(fname):
        logger.error("Inputfile %s does not exist!", fname)
        sys.exit(1)

    for i, line in enumerate(open(fname)):
        if line.startswith("#") or line.startswith("track") or line.startswith("browser"):
            # comment or BED specific stuff
            pass
        else:
            vals = line.strip().split("\t")
            if len(vals) < 3:
                logger.error("Expecting tab-seperated values (chromosome<tab>start<tab>end) on line %s of file %s", i + 1, fname)
                sys.exit(1)
            try:
                start, end = int(vals[1]), int(vals[2])
            except ValueError:
                logger.error("No valid integer coordinates on line %s of file %s", i + 1, fname)
                sys.exit(1)
            if len(vals) > 3:
                try:
                    float(vals[3])
                except ValueError:
                    pass
                    #self.logger.warn("No numerical value in column 4 on line %s of file %s, ignoring..." % (i + 1, file))


def check_denovo_input(inputfile, params):

    genome = params["genome"]
    background = params["background"]
    
    input_type = "BED"
    # If we can load it as fasta then it is a fasta, yeh?
    try:
        Fasta(inputfile)
        logger.debug("Inputfile is a FASTA file")
        input_type = "FASTA"
    except Exception:
        # Leave it to BED
        pass

    if input_type == "FASTA":
        valid_bg = FA_VALID_BGS    
    elif input_type == "BED":
        valid_bg = BED_VALID_BGS    
        if "genomic" in background:
            check_genome(genome)
        # is it a valid bed-file etc.
        check_bed_file(inputfile)    # bed-specific
    
    for bg in background:
        if not bg in valid_bg:
            logger.info("Input type is %s, ignoring background type '%s'", 
                            input_type, bg)
        background = [bg for bg in background if bg in valid_bg]

    if len(background) == 0:
        logger.error("No valid backgrounds specified!")
        sys.exit(1)

    return input_type, background


