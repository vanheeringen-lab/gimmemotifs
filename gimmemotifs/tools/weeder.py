from .motifprogram import MotifProgram
import os

import shutil
from subprocess import Popen, PIPE
from tempfile import NamedTemporaryFile

from gimmemotifs.motif import read_motifs


class Weeder(MotifProgram):

    """
    Predict motifs using Weeder.

    Reference:

    """

    def __init__(self):
        self.name = "Weeder"
        self.cmd = "weeder2"
        self.use_width = False
        self.default_params = {"organism": "hg19", "single": False}

    def _run_program(self, bin, fastafile, params=None):
        """
        Run Weeder and predict motifs from a FASTA file.

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
        params = self._parse_params(params)

        organism = params["organism"]
        weeder_organisms = {
            "hg18": "HS",
            "hg19": "HS",
            "hg38": "HS",
            "mm9": "MM",
            "mm10": "MM",
            "dm3": "DM",
            "dm5": "DM",
            "dm6": "DM",
            "yeast": "SC",
            "sacCer2": "SC",
            "sacCer3": "SC",
            "TAIR10": "AT",
            "TAIR11": "AT",
        }
        weeder_organism = weeder_organisms.get(organism, "HS")

        tmp = NamedTemporaryFile(dir=self.tmpdir)
        name = tmp.name
        tmp.close()
        shutil.copy(fastafile, name)
        fastafile = name

        cmd = "{} -f {} -O {}".format(self.cmd, fastafile, weeder_organism)

        if params["single"]:
            cmd += " -ss"

        # print cmd
        stdout, stderr = "", ""
        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE, cwd=self.tmpdir)
        out, err = p.communicate()
        stdout += out.decode()
        stderr += err.decode()

        motifs = []
        if os.path.exists(fastafile + ".matrix.w2"):
            f = open(fastafile + ".matrix.w2")
            motifs = self.parse(f)
            f.close()

        for m in motifs:
            m.id = "{}_{}".format(self.name, m.id.split("\t")[0])

        for ext in [".w2", ".matrix.w2"]:
            if os.path.exists(fastafile + ext):
                os.unlink(fastafile + ext)

        return motifs, stdout, stderr

    def parse(self, fo):
        """
        Convert Weeder output to motifs

        Parameters
        ----------
        fo : file-like
            File object containing Weeder output.

        Returns
        -------
        motifs : list
            List of Motif instances.
        """
        return read_motifs(fo, fmt="jaspar")
