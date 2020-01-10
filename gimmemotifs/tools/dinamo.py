from .motifprogram import MotifProgram
import os
from subprocess import Popen, PIPE

from gimmemotifs.motif import read_motifs


class Dinamo(MotifProgram):

    """
    Predict motifs using DiNAMO.

    Reference: Saad et al., 2018, doi: 10.1186/s12859-018-2215-1
    """

    def __init__(self):
        self.name = "DiNAMO"
        self.cmd = "dinamo"
        self.use_width = True
        self.default_params = {"background": None, "number": 10, "width": 10}

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
            raise ValueError("DiNAMO needs a background file.")

        return prm

    def _run_program(self, bin, fastafile, params=None):
        """
        Run DiNAMO and predict motifs from a FASTA file.

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

        outfile = os.path.join(self.tmpdir, "motifs.meme")
        stdout = ""
        stderr = ""

        cmd = "%s -pf %s -nf %s -l %s -o %s" % (
            bin,
            fastafile,
            params["background"],
            params["width"],
            outfile,
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
