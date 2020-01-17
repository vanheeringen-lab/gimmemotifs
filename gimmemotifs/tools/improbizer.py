from .motifprogram import MotifProgram
import os
import re
from subprocess import Popen, PIPE

from gimmemotifs.motif import Motif


class Improbizer(MotifProgram):

    """
    Predict motifs using Improbizer.

    Reference:
    """

    def __init__(self):
        self.name = "Improbizer"
        self.cmd = "ameme"
        self.use_width = False
        self.default_params = {"background": None, "number": 10}

    def _parse_params(self, params=None):
        """
        Parse parameters.

        Combine default and user-defined parameters.
        """
        prm = super()._parse_params(params, needs_background=True)
        prm["outfile"] = os.path.join(self.tmpdir, "improbizer.out.html")
        return prm

    def _run_program(self, bin, fastafile, params=None):
        """
        Run Improbizer and predict motifs from a FASTA file.

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

        current_path = os.getcwd()
        os.chdir(self.tmpdir)

        stdout = ""
        stderr = ""
        cmd = "%s good=%s bad=%s numMotifs=%s > %s" % (
            bin,
            fastafile,
            params["background"],
            params["number"],
            params["outfile"],
        )
        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()
        stdout += out.decode()
        stderr += err.decode()

        motifs = []
        if os.path.exists(params["outfile"]):
            f = open(params["outfile"])
            motifs = self.parse(f)
            f.close()

        os.chdir(current_path)

        return motifs, stdout, stderr

    def parse(self, fo):
        """
        Convert Improbizer output to motifs

        Parameters
        ----------
        fo : file-like
            File object containing Improbizer output.

        Returns
        -------
        motifs : list
            List of Motif instances.
        """
        motifs = []
        p = re.compile(r"\d+\s+@\s+\d+\.\d+\s+sd\s+\d+\.\d+\s+(\w+)$")

        line = fo.readline()
        while line and line.find("Color") == -1:
            m = p.search(line)
            if m:
                pwm_data = {}
                for _i in range(4):
                    vals = [x.strip() for x in fo.readline().strip().split(" ") if x]
                    pwm_data[vals[0].upper()] = vals[1:]
                pwm = []
                for i in range(len(pwm_data["A"])):
                    pwm.append([float(pwm_data[x][i]) for x in ["A", "C", "G", "T"]])
                motifs.append(Motif(pwm))
                motifs[-1].id = "%s_%s" % (self.name, m.group(1))
            line = fo.readline()

        return motifs
