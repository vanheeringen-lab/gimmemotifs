# Copyright (c) 2009-2019 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.

""" Odds and ends that for which I didn't (yet) find another place """
# Python imports
import os
import re
import sys
import hashlib
import logging
import mmap
import random
import tempfile
import requests
from io import TextIOWrapper
from functools import singledispatch
from subprocess import Popen
from tempfile import NamedTemporaryFile
from shutil import copyfile

# External imports
import pyfaidx
from scipy import special
import numpy as np
import pybedtools

from genomepy import Genome
from Bio.SeqIO.FastaIO import SimpleFastaParser

# gimme imports
from gimmemotifs.fasta import Fasta
from gimmemotifs.plot import plot_histogram
from gimmemotifs.rocmetrics import ks_pvalue
from gimmemotifs.config import MotifConfig


logger = logging.getLogger("gimme.utils")

# pylint: disable=no-member
lgam = special.gammaln


def rc(seq):
    """ Return reverse complement of sequence """
    d = str.maketrans("actgACTG", "tgacTGAC")
    return seq[::-1].translate(d)


def narrowpeak_to_bed(inputfile, bedfile, size=0):
    """Convert narrowPeak file to BED file."""
    p = re.compile(r"^(#|track|browser)")
    warn_no_summit = True
    with open(bedfile, "w") as f_out:
        with open(inputfile) as f_in:
            for line in f_in:
                if p.search(line):
                    continue
                vals = line.strip().split("\t")
                start, end = int(vals[1]), int(vals[2])

                if size > 0:
                    summit = int(vals[9])
                    if summit == -1:
                        if warn_no_summit:
                            logger.warn(
                                "No summit present in narrowPeak file, "
                                "using the peak center."
                            )
                            warn_no_summit = False
                        summit = (end - start) // 2

                    start = start + summit - (size // 2)
                    end = start + size
                f_out.write("{}\t{}\t{}\t{}\n".format(vals[0], start, end, vals[6]))


def pfmfile_location(infile):
    config = MotifConfig()

    if infile is None:
        infile = config.get_default_params().get("motif_db", None)
        if infile is None:
            raise ValueError(
                "No motif file was given and no default "
                "database specified in the config file."
            )

    if isinstance(infile, str):
        if not os.path.exists(infile):
            motif_dir = config.get_motif_dir()
            checkfile = os.path.join(motif_dir, infile)
            if os.path.exists(checkfile):
                infile = checkfile
            else:
                for ext in [".pfm", ".pwm"]:
                    if os.path.exists(checkfile + ext):
                        infile = checkfile + ext
                        break
            if not os.path.exists(infile):
                raise ValueError("Motif file {} not found".format(infile))

    return infile


def get_jaspar_motif_info(motif_id):
    query_url = "http://jaspar.genereg.net/api/v1/matrix/{}?format=json"
    result = requests.get(query_url.format(motif_id))

    if not result.ok:
        result.raise_for_status()
        sys.exit()

    return result.json()


def phyper_single(k, good, bad, N):

    return np.exp(
        lgam(good + 1)
        - lgam(good - k + 1)
        - lgam(k + 1)
        + lgam(bad + 1)
        - lgam(bad - N + k + 1)
        - lgam(N - k + 1)
        - lgam(bad + good + 1)
        + lgam(bad + good - N + 1)
        + lgam(N + 1)
    )


def phyper(k, good, bad, N):
    """Current hypergeometric implementation in scipy is broken,
    so here's the correct version.
    """
    pvalues = [phyper_single(x, good, bad, N) for x in range(k + 1, N + 1)]
    return np.sum(pvalues)


def divide_file(fname, sample, rest, fraction, abs_max):
    with open(fname) as f:
        lines = f.readlines()
    # random.seed()
    random.shuffle(lines)

    x = int(fraction * len(lines))
    if x > abs_max:
        x = abs_max

    tmp = tempfile.NamedTemporaryFile(mode="w", delete=False)

    # Fraction as sample
    for line in lines[:x]:
        tmp.write(line)
    tmp.flush()

    # Make sure it is sorted for tools that use this information (MDmodule)
    stdout, stderr = Popen(
        "sort -k4gr %s > %s" % (tmp.name, sample), shell=True
    ).communicate()

    tmp.close()

    if stderr:
        print("Something went wrong.\nstdout: {}\nstderr; {}".format(stdout, stderr))
        sys.exit()

    # Rest
    f = open(rest, "w")
    for line in lines[x:]:
        f.write(line)
    f.close()

    # if os.path.exists(tmp.name):
    #    os.unlink(tmp.name)
    return x, len(lines[x:])


def divide_fa_file(fname, sample, rest, fraction, abs_max):
    fa = Fasta(fname)
    ids = fa.ids[:]

    x = int(fraction * len(ids))
    if x > abs_max:
        x = abs_max

    sample_seqs = random.sample(ids, x)

    # Rest
    f_sample = open(sample, "w")
    f_rest = open(rest, "w")
    for name, seq in fa.items():
        if name in sample_seqs:
            f_sample.write(">%s\n%s\n" % (name, seq))
        else:
            f_rest.write(">%s\n%s\n" % (name, seq))
    f_sample.close()
    f_rest.close()

    return x, len(ids[x:])


def write_equalsize_bedfile(bedfile, size, outfile):
    """Read input from <bedfile>, set the size of all entries to <size> and
    write the result to <outfile>.
    Input file needs to be in BED or WIG format."""
    if size is None or size <= 0:
        copyfile(bedfile, outfile)
        return

    BUFSIZE = 10000
    f = open(bedfile)
    out = open(outfile, "w")
    lines = f.readlines(BUFSIZE)
    line_count = 0
    while lines:
        for line in lines:
            line_count += 1
            if (
                not line.startswith("#")
                and not line.startswith("track")
                and not line.startswith("browser")
            ):
                vals = line.strip().split("\t")
                try:
                    start, end = int(vals[1]), int(vals[2])
                except ValueError:
                    print(
                        "Error on line %s while reading %s. "
                        "Is the file in BED or WIG format?" % (line_count, bedfile)
                    )
                    sys.exit(1)

                start = (start + end) // 2 - (size // 2)
                # This shifts the center, but ensures the size is identical...
                # maybe not ideal
                if start < 0:
                    start = 0
                end = start + size
                # Keep all the other information in the bedfile if it's there
                if len(vals) > 3:
                    out.write(
                        "%s\t%s\t%s\t%s\n" % (vals[0], start, end, "\t".join(vals[3:]))
                    )
                else:
                    out.write("%s\t%s\t%s\n" % (vals[0], start, end))
        lines = f.readlines(BUFSIZE)

    out.close()
    f.close()


def median_bed_len(bedfile):
    f = open(bedfile)
    lengths = []
    for i, line in enumerate(f.readlines()):
        if not (line.startswith("browser") or line.startswith("track")):
            vals = line.split("\t")
            try:
                lengths.append(int(vals[2]) - int(vals[1]))
            except ValueError:
                sys.stderr.write(
                    "Error in line %s: "
                    "coordinates in column 2 and 3 need to be integers!\n" % (i)
                )
                sys.exit(1)
    f.close()
    return np.median(lengths)


def motif_localization(fastafile, motif, size, outfile, cutoff=0.9):
    NR_HIST_MATCHES = 100

    matches = motif.pwm_scan(Fasta(fastafile), cutoff=cutoff, nreport=NR_HIST_MATCHES)
    if len(matches) > 0:
        ar = []
        for a in matches.values():
            ar += a
        matches = np.array(ar)
        p = ks_pvalue(matches, size - len(motif))
        plot_histogram(
            matches - size / 2 + len(motif) / 2,
            outfile,
            xrange=(-size / 2, size / 2),
            breaks=21,
            title="%s (p=%0.2e)" % (motif.id, p),
            xlabel="Position",
        )
        return motif.id, p
    else:
        return motif.id, 1.0


def parse_cutoff(motifs, cutoff, default=0.9):
    """Provide either a file with one cutoff per motif or a single cutoff
    returns a hash with motif id as key and cutoff as value
    """

    cutoffs = {}
    if os.path.isfile(str(cutoff)):
        for i, line in enumerate(open(cutoff)):
            if line != "Motif\tScore\tCutoff\n":
                try:
                    motif, _, c = line.strip().split("\t")
                    c = float(c)
                    cutoffs[motif] = c
                except Exception as e:
                    sys.stderr.write(
                        "Error parsing cutoff file, line {0}: {1}\n".format(e, i + 1)
                    )
                    sys.exit(1)
    else:
        for motif in motifs:
            cutoffs[motif.id] = float(cutoff)

    for motif in motifs:
        if motif.id not in cutoffs:
            sys.stderr.write(
                "No cutoff found for {0}, using default {1}\n".format(motif.id, default)
            )
            cutoffs[motif.id] = default
    return cutoffs


def _treesort(order, nodeorder, nodecounts, tree):
    # From the Pycluster library, Michiel de Hoon
    # Find the order of the nodes consistent with the hierarchical clustering
    # tree, taking into account the preferred order of nodes.
    nNodes = len(tree)
    nElements = nNodes + 1
    neworder = np.zeros(nElements)
    clusterids = np.arange(nElements)
    for i in range(nNodes):
        i1 = tree[i].left
        i2 = tree[i].right
        if i1 < 0:
            order1 = nodeorder[-i1 - 1]
            count1 = nodecounts[-i1 - 1]
        else:
            order1 = order[i1]
            count1 = 1
        if i2 < 0:
            order2 = nodeorder[-i2 - 1]
            count2 = nodecounts[-i2 - 1]
        else:
            order2 = order[i2]
            count2 = 1
        # If order1 and order2 are equal, their order is determined
        # by the order in which they were clustered
        if i1 < i2:
            if order1 < order2:
                increase = count1
            else:
                increase = count2
            for j in range(nElements):
                clusterid = clusterids[j]
                if clusterid == i1 and order1 >= order2:
                    neworder[j] += increase
                if clusterid == i2 and order1 < order2:
                    neworder[j] += increase
                if clusterid == i1 or clusterid == i2:
                    clusterids[j] = -i - 1
        else:
            if order1 <= order2:
                increase = count1
            else:
                increase = count2
            for j in range(nElements):
                clusterid = clusterids[j]
                if clusterid == i1 and order1 > order2:
                    neworder[j] += increase
                if clusterid == i2 and order1 <= order2:
                    neworder[j] += increase
                if clusterid == i1 or clusterid == i2:
                    clusterids[j] = -i - 1
    return np.argsort(neworder)


def number_of_seqs_in_file(fname):
    try:
        fa = Fasta(fname)
        return len(fa)
    except Exception:
        pass

    try:
        bed = pybedtools.BedTool(fname)
        return len([x for x in bed])
    except Exception:
        pass

    sys.stderr.write("unknown filetype {}\n".format(fname))
    sys.exit(1)


def determine_file_type(fname):
    """
    Detect file type.

    The following file types are supported:
    BED, narrowPeak, FASTA, list of chr:start-end regions
    If the extension is bed, fa, fasta or narrowPeak, we will believe this
    without checking!

    Parameters
    ----------
    fname : str
        File name.

    Returns
    -------
    filetype : str
        Filename in lower-case.
    """
    if not (isinstance(fname, str)):
        raise ValueError(f"{fname} is not a string, does not represent a file name")

    if not os.path.exists(fname):
        raise ValueError(f"File {fname} does not exist!")

    if not os.path.isfile(fname):
        raise ValueError(f"{fname} is not a file!")

    ext = os.path.splitext(fname)[1].lower()
    if ext in [".bed"]:
        return "bed"
    elif ext in [".fa", ".fasta"]:
        return "fasta"
    elif ext in [".narrowpeak"]:
        return "narrowpeak"

    try:
        Fasta(fname)
        return "fasta"
    except Exception:
        pass
    # Read first line that is not a comment or an UCSC-specific line
    p = re.compile(r"^(#|track|browser)")
    with open(fname) as f:
        for line in f.readlines():
            line = line.strip()
            if not p.search(line):
                break
    region_p = re.compile(r"^(.+):(\d+)-(\d+)$")
    if region_p.search(line):
        return "region"
    else:
        vals = line.split("\t")
        if len(vals) >= 3:
            try:
                _, _ = int(vals[1]), int(vals[2])
            except ValueError:
                return "unknown"

            if len(vals) == 10:
                try:
                    _, _ = int(vals[4]), int(vals[9])
                    return "narrowpeak"
                except ValueError:
                    # As far as I know there is no 10-column BED format
                    return "unknown"
            return "bed"

    # Catch-all
    return "unknown"


def get_seqs_type(seqs):
    """
    automagically determine input type
    the following types are detected:
        - Fasta object
        - FASTA file
        - list of regions
        - region file
        - BED file
    """
    region_p = re.compile(r"^([^\s:]+\@)?(.+):(\d+)-(\d+)$")
    if isinstance(seqs, Fasta):
        return "fasta"
    elif isinstance(seqs, list) or isinstance(seqs, np.ndarray):
        if len(seqs) == 0:
            raise ValueError("empty list of sequences to scan")
        else:
            if region_p.search(seqs[0]):
                return "regions"
            else:
                raise ValueError("unknown region type")
    elif isinstance(seqs, str):
        if os.path.isfile(seqs):
            ftype = determine_file_type(seqs)
            if ftype == "unknown":
                raise ValueError("unknown type")
            elif ftype == "narrowpeak":
                raise ValueError("narrowPeak not yet supported in this function")
            else:
                return ftype + "file"
        else:
            raise ValueError("no file found with name {}".format(seqs))
    else:
        raise ValueError("unknown type {}".format(type(seqs).__name__))


# Regular expression to check for region (chr:start-end or genome@chr:start-end)
region_p = re.compile(r"^[^@]+@([^\s]+):(\d+)-(\d+)$")


def _check_minsize(fa, minsize):
    """
    Raise ValueError if there is any sequence that is shorter than minsize.
    If minsize is None the size will not be checked.
    """
    if minsize is None:
        return fa

    for name, seq in fa.items():
        if len(seq) < minsize:
            raise ValueError(f"sequence {name} is shorter than {minsize}")

    return fa


def _genomepy_convert(to_convert, genome, minsize=None):
    """
    Convert a variety of inputs using track2fasta().
    """
    if genome is None:
        raise ValueError("input file is not a FASTA file, need a genome!")

    if isinstance(genome, Genome):
        g = genome
    else:
        g = Genome(genome)

    tmpfile = NamedTemporaryFile()
    g.track2fasta(to_convert, tmpfile.name)

    fa = as_seqdict(tmpfile.name)
    return _check_minsize(fa, minsize)


def _as_seqdict_genome_regions(regions, minsize=None):
    """
    Accepts list of regions where the genome is encoded in the region,
    using the genome@chrom:start-end format.
    """
    genomic_regions = {}
    for region in regions:
        genome, region = region.split("@")
        if genome not in genomic_regions:
            Genome(genome)
            genomic_regions[genome] = []
        genomic_regions[genome].append(region)

    tmpfa = NamedTemporaryFile(mode="w", delete=False)
    for genome, g_regions in genomic_regions.items():
        g = Genome(genome)

        fa = g.track2fasta(g_regions)

        for seq in fa:
            seq.name = f"{genome}@{seq.name}"
            print(seq.__repr__(), file=tmpfa)

    tmpfa.flush()

    # Open tempfile and restore original sequence order
    fa = as_seqdict(tmpfa.name)
    fa = {region: fa[region] for region in regions}
    return _check_minsize(fa, minsize)


@singledispatch
def as_seqdict(to_convert, genome=None, minsize=None):
    """
    Convert input to a dictionary with name as key and sequence as value.

    If the input contains genomic coordinates, the genome needs to be
    specified. If minsize is specified all sequences will be checked if they
    are not shorter than minsize. If regions (or a region file) are used as
    the input, the genome can optionally be specified in the region using the
    following format: genome@chrom:start-end.

    Current supported input types include:
    * FASTA, BED and region files.
    * List or numpy.ndarray of regions.
    * pyfaidx.Fasta object.
    * pybedtools.BedTool object.

    Parameters
    ----------
    to_convert : list, str, pyfaidx.Fasta or pybedtools.BedTool
        Input to convert to FASTA-like dictionary

    genome : str, optional
        Genomepy genome name.

    minsize : int or None, optional
        If specified, check if all sequences have at least size minsize.

    Returns
    -------
        dict with sequence names as key and sequences as value.
    """
    raise NotImplementedError(f"Not implement for {type(to_convert)}")


@as_seqdict.register(list)
def _as_seqdict_list(to_convert, genome=None, minsize=None):
    """
    Accepts list of regions as input.
    """
    if region_p.match(to_convert[0]):
        return _as_seqdict_genome_regions(to_convert, minsize)

    return _genomepy_convert(to_convert, genome, minsize)


@as_seqdict.register(TextIOWrapper)
def _as_seqdict_file_object(to_convert, genome=None, minsize=None):
    """
    Accepts file object as input, should be a FASTA file.
    """
    fa = {x: y for x, y in SimpleFastaParser(to_convert)}
    return _check_minsize(fa, minsize)


@as_seqdict.register(str)
def _as_seqdict_filename(to_convert, genome=None, minsize=None):
    """
    Accepts filename as input.
    """
    if not os.path.exists(to_convert):
        raise ValueError("Assuming filename, but it does not exist")

    f = open(to_convert)
    fa = as_seqdict(f)

    if any(fa):
        return _check_minsize(fa, minsize)

    with open(to_convert) as f:
        line = ""
        while True:
            line = f.readline()
            if line == "":
                break
            if not line.startswith("#"):
                break

        if line == "":
            raise IOError(f"empty file {to_convert}")

        if region_p.match(line.strip()):
            regions = [myline.strip() for myline in [line] + f.readlines()]
            return _as_seqdict_genome_regions(regions, minsize=None)

    # Biopython parser resulted in empty dict
    # Assuming it's a BED or region file
    return _genomepy_convert(to_convert, genome, minsize)


@as_seqdict.register(pyfaidx.Fasta)
def _as_seqdict_pyfaidx(to_convert, genome=None, minsize=None):
    """
    Accepts pyfaidx.Fasta object as input.
    """
    fa = {k: str(v) for k, v in to_convert.items()}
    return _check_minsize(fa, minsize)


@as_seqdict.register(pybedtools.BedTool)
def _as_seqdict_bedtool(to_convert, genome=None, minsize=None):
    """
    Accepts pybedtools.BedTool as input.
    """
    return _genomepy_convert(
        ["{}:{}-{}".format(*f[:3]) for f in to_convert], genome, minsize
    )


@as_seqdict.register(np.ndarray)
def _as_seqdict_array(to_convert, genome=None, minsize=None):
    """
    Accepts numpy.ndarray with regions as input.
    """
    return as_seqdict(list(to_convert), genome, minsize)


def as_fasta(to_convert, genome=None, minsize=None):
    if isinstance(to_convert, Fasta):
        return to_convert

    return Fasta(fdict=as_seqdict(to_convert, genome, minsize))


def file_checksum(fname):
    """Return md5 checksum of file.

    Note: only works for files < 4GB.

    Parameters
    ----------
    filename : str
        File used to calculate checksum.

    Returns
    -------
        checkum : str
    """
    size = os.path.getsize(fname)
    with open(fname, "r+") as f:
        checksum = hashlib.md5(mmap.mmap(f.fileno(), size)).hexdigest()
    return checksum


def join_max(a, length, sep="", suffix=""):
    lengths = [len(x) for x in a]
    total = 0
    for i, size in enumerate(lengths + [0]):
        if total > (length - len(suffix)):
            return sep.join(a[: i - 1]) + suffix
        if i > 0:
            total += 1
        total += size
    return sep.join(a)


def check_genome(genome):
    """Check if genome is a valid FASTA file or genomepy genome genome.

    Parameters
    ----------
    genome : str
        Genome name or file to check.

    Returns
    -------
    is_genome : bool
    """
    try:
        Genome(genome)
        return True
    except Exception:
        pass
    return False
