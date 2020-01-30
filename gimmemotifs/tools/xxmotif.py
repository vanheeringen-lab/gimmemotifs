from .motifprogram import MotifProgram
import os
from subprocess import Popen, PIPE


class XXmotif(MotifProgram):

    """
    Predict motifs using XXmotif.

    Reference:
    """

    def __init__(self):
        self.name = "XXmotif"
        self.cmd = "XXmotif"
        self.use_width = False
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
        prm = super()._parse_params(params)

        if prm["background"]:
            # Absolute path, just to be sure
            prm["background"] = " --negSet {0} ".format(prm["background"])

        prm["strand"] = ""
        if not prm["single"]:
            prm["strand"] = " --revcomp "

        return prm

    def _run_program(self, bin, fastafile, params=None):
        """
        Run XXmotif and predict motifs from a FASTA file.

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

        outfile = os.path.join(
            self.tmpdir, os.path.basename(fastafile.replace(".fa", ".pwm"))
        )

        stdout = ""
        stderr = ""

        cmd = "%s %s %s --localization --batch %s %s" % (
            bin,
            self.tmpdir,
            fastafile,
            params["background"],
            params["strand"],
        )

        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()
        stdout += out.decode()
        stderr += err.decode()

        motifs, stdout, stderr = self._read_and_label_motifs(
            outfile, stdout, stderr, fmt="xxmotif"
        )

        return motifs, stdout, stderr
