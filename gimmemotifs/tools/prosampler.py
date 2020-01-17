from .motifprogram import MotifProgram
import os
from subprocess import Popen, PIPE


class ProSampler(MotifProgram):

    """
    Predict motifs using ProSampler.

    Reference: Li et al., 2019, doi: 10.1093/bioinformatics/btz290
    """

    def __init__(self):
        self.name = "ProSampler"
        self.cmd = "ProSampler"
        self.use_width = False
        self.default_params = {"single": False, "background": None}

    def _parse_params(self, params=None):
        """
        Parse parameters.

        Combine default and user-defined parameters.
        """
        prm = super()._parse_params(params, needs_background=True)
        prm["strand"] = " -p 2 "
        if prm["single"]:
            prm["strand"] = " -p 1 "

        return prm

    def _run_program(self, bin, fastafile, params=None):
        """
        Run ProSampler and predict motifs from a FASTA file.

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

        outfile = os.path.join(self.tmpdir, "ProSampler.meme")
        stdout = ""
        stderr = ""

        cmd = "%s -i %s -b %s -o %s %s" % (
            bin,
            fastafile,
            params["background"],
            os.path.join(self.tmpdir, "ProSampler"),
            params["strand"],
        )

        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE, cwd=self.tmpdir)
        out, err = p.communicate()
        stdout += out.decode()
        stderr += err.decode()

        motifs, stdout, stderr = self._read_and_label_motifs(
            outfile, stdout, stderr, fmt="meme"
        )

        return motifs, stdout, stderr
