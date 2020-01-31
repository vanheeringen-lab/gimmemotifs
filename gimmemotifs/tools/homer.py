from .motifprogram import MotifProgram
import os
from subprocess import Popen, PIPE
from tempfile import NamedTemporaryFile

from gimmemotifs.motif import read_motifs


class Homer(MotifProgram):

    """
    Predict motifs using Homer.

    Reference: Heinz et al, 2010; PMID: 20513432
    """

    def __init__(self):
        self.name = "Homer"
        self.cmd = "homer2"
        self.use_width = True
        self.default_params = {
            "single": False,
            "background": None,
            "analysis": "medium",
            "number": 5,
            "width": 10,
        }

    def _parse_params(self, params=None):
        """
        Parse parameters.

        Combine default and user-defined parameters.
        """
        prm = super()._parse_params(params, needs_background=True)
        prm["strand"] = ""
        if prm["single"]:
            prm["strand"] = " -strand + "

        return prm

    def _run_program(self, bin, fastafile, params=None):
        """
        Run Homer and predict motifs from a FASTA file.

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

        outfile = NamedTemporaryFile(
            mode="w", dir=self.tmpdir, prefix="homer_w{}.".format(params["width"])
        ).name

        cmd = "%s denovo -i %s -b %s -len %s -S %s %s -o %s -p 8" % (
            bin,
            fastafile,
            params["background"],
            params["width"],
            params["number"],
            params["strand"],
            outfile,
        )

        stderr = ""
        stdout = "Running command:\n{}\n".format(cmd)

        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE, cwd=self.tmpdir)
        out, err = p.communicate()
        stdout += out.decode()
        stderr += err.decode()

        motifs = []

        if os.path.exists(outfile):
            motifs = read_motifs(outfile, fmt="pwm")
            for i, m in enumerate(motifs):
                m.id = "{}_{}_{}".format(self.name, params["width"], i + 1)

        return motifs, stdout, stderr
