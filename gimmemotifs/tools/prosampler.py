from .motifprogram import MotifProgram
import os
from subprocess import Popen, PIPE

from gimmemotifs.motif import read_motifs


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
        prm = self.default_params.copy()
        if params is not None:
            prm.update(params)

        if prm["background"]:
            # Absolute path, just to be sure
            prm["background"] = os.path.abspath(prm["background"])
        else:
            raise ValueError("ProSampler needs a background file.")

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

        motifs = []

        if os.path.exists(outfile):
            motifs = read_motifs(outfile, fmt="meme")
            for m in motifs:
                m.id = "{0}_{1}".format(self.name, m.id)
        else:
            stdout += "\nMotif file {0} not found!\n".format(outfile)
            stderr += "\nMotif file {0} not found!\n".format(outfile)

        return motifs, stdout, stderr
