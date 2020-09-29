from .motifprogram import MotifProgram
import io
import os
import re
from subprocess import Popen, PIPE
from tempfile import NamedTemporaryFile

from gimmemotifs.motif import Motif


class Meme(MotifProgram):

    """
    Predict motifs using MEME.

    Reference:
    """

    def __init__(self):
        self.name = "MEME"
        self.cmd = "meme"
        self.use_width = True

    def _run_program(self, bin, fastafile, params=None):
        """
        Run MEME and predict motifs from a FASTA file.

        Parameters
        ----------
        bin : str
            Command used to run the tool.

        fastafile : str
            Name of the FASTA input file.

        params : dict, optional
            Optional parameters. For some of the tools required parameters
            are passed using this dictionary.

        Returns
        -------
        motifs : list of Motif instances
            The predicted motifs.

        stdout : str
            Standard out of the tool.

        stderr : str
            Standard error of the tool.
        """
        default_params = {"width": 10, "single": False, "number": 10}
        if params is not None:
            default_params.update(params)

        tmp = NamedTemporaryFile(dir=self.tmpdir)

        strand = "-revcomp"
        width = default_params["width"]
        number = default_params["number"]

        cmd = [
            bin,
            fastafile,
            "-text",
            "-dna",
            "-nostatus",
            "-mod",
            "zoops",
            "-nmotifs",
            "%s" % number,
            "-w",
            "%s" % width,
            "-maxsize",
            "10000000",
        ]
        if not default_params["single"]:
            cmd.append(strand)

        # Fix to run in Docker
        env = os.environ.copy()
        env["OMPI_MCA_plm_rsh_agent"] = "sh"

        p = Popen(cmd, bufsize=1, stderr=PIPE, stdout=PIPE, env=env)
        stdout, stderr = p.communicate()

        motifs = []
        motifs = self.parse(io.StringIO(stdout.decode()))

        # Delete temporary files
        tmp.close()

        return motifs, stdout, stderr

    def parse(self, fo):
        """
        Convert MEME output to motifs

        Parameters
        ----------
        fo : file-like
            File object containing MEME output.

        Returns
        -------
        motifs : list
            List of Motif instances.
        """
        motifs = []
        nucs = {"A": 0, "C": 1, "G": 2, "T": 3}

        p = re.compile(r"MOTIF.+MEME-(\d+)\s*width\s*=\s*(\d+)\s+sites\s*=\s*(\d+)")
        pa = re.compile(r"\)\s+([A-Z]+)")
        line = fo.readline()
        while line:
            m = p.search(line)
            align = []
            pfm = None
            if m:
                # print(m.group(0))
                id = "%s_%s_w%s" % (self.name, m.group(1), m.group(2))
                while not line.startswith("//"):
                    ma = pa.search(line)
                    if ma:
                        # print(ma.group(0))
                        match = ma.group(1)
                        align.append(match)
                        if not pfm:
                            pfm = [[0 for x in range(4)] for x in range(len(match))]
                        for pos in range(len(match)):
                            if match[pos] in nucs:
                                pfm[pos][nucs[match[pos]]] += 1
                            else:
                                for i in range(4):
                                    pfm[pos][i] += 0.25

                    line = fo.readline()

                motifs.append(Motif(pfm[:]))
                motifs[-1].id = id
                motifs[-1].align = align[:]
            line = fo.readline()

        return motifs
