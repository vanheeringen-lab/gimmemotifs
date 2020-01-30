from .motifprogram import MotifProgram
import os
import re
from subprocess import Popen, PIPE

from gimmemotifs.motif import Motif


class BioProspector(MotifProgram):

    """
    Predict motifs using BioProspector.

    Reference:
    """

    def __init__(self):
        self.name = "BioProspector"
        self.cmd = "BioProspector"
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
        prm["strand"] = 2
        if prm["single"]:
            prm["strand"] = 1

        return prm

    def _run_program(self, bin, fastafile, params=None):
        """
        Run BioProspector and predict motifs from a FASTA file.

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

        outfile = os.path.join(self.tmpdir, "bioprospector.out")

        stdout = ""
        stderr = ""

        cmd = "%s -i %s -W %s -d %s -b %s -r %s -o %s" % (
            bin,
            fastafile,
            params["width"],
            params["strand"],
            params["background"],
            params["number"],
            outfile,
        )

        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()
        stdout += out.decode()
        stderr += err.decode()

        motifs = []

        if os.path.exists(outfile):
            with open(outfile) as f:
                motifs = self.parse(f)

        return motifs, stdout, stderr

    def parse(self, fo):
        """
        Convert BioProspector output to motifs

        Parameters
        ----------
        fo : file-like
            File object containing BioProspector output.

        Returns
        -------
        motifs : list
            List of Motif instances.
        """
        motifs = []

        p = re.compile(r"^\d+\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)")
        pwm = []
        motif_id = ""
        for line in fo.readlines():
            if line.startswith("Motif #"):
                if pwm:
                    m = Motif(pwm)
                    m.id = "BioProspector_w%s_%s" % (len(m), motif_id)
                    motifs.append(m)
                motif_id = line.split("#")[1].split(":")[0]
                pwm = []
            else:
                m = p.search(line)
                if m:
                    pwm.append([float(m.group(x)) / 100.0 for x in range(1, 5)])

        if pwm:
            m = Motif(pwm)
            m.id = "BioProspector_w%s_%s" % (len(m), motif_id)
            motifs.append(m)
        return motifs
