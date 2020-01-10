from .motifprogram import MotifProgram
import os
from subprocess import Popen, PIPE

from gimmemotifs.motif import read_motifs


class Yamda(MotifProgram):

    """
    Predict motifs using YAMDA.

    Reference: Quang et al., 2018, https://doi.org/10.1093/bioinformatics/bty396
    """

    def __init__(self):
        self.name = "YAMDA"
        self.cmd = "run_em.py"
        self.use_width = True
        self.default_params = {
            "single": False,
            "background": None,
            "number": 10,
            "width": 20,
        }

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
            raise ValueError("YAMDA needs a background file.")

        prm["strand"] = ""
        if not prm["single"]:
            prm["strand"] = " --revcomp "

        return prm

    def _run_program(self, bin, fastafile, params=None):
        """
        Run YAMDA and predict motifs from a FASTA file.

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

        outfile = os.path.join(self.tmpdir, "motifs.txt")
        stdout = ""
        stderr = ""

        cmd = "%s -i %s -j %s -n %s -w %s -oc %s  %s" % (
            bin,
            fastafile,
            params["background"],
            params["number"],
            params["width"],
            self.tmpdir,
            params["strand"],
        )

        print(cmd)
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
