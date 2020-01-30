from .motifprogram import MotifProgram
import os
from subprocess import Popen, PIPE

from gimmemotifs.motif import Motif


class Rpmcmc(MotifProgram):

    """
    Predict motifs using RPMCMC.

    Reference: Ikebata & Yoshida, 2015, 10.1093/bioinformatics/btv017
    """

    def __init__(self):
        self.name = "RPMCMC"
        self.cmd = "multi_motif_finder"
        self.use_width = False
        self.default_params = {}

    def _run_program(self, bin, fastafile, params=None):
        """
        Run RPMCMC and predict motifs from a FASTA file.

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

        outfile = os.path.join(self.tmpdir, "pfm")
        stdout = ""
        stderr = ""

        cmd = "ulimit -s unlimited && %s -d %s -od ./" % (bin, fastafile)

        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE, cwd=self.tmpdir)
        out, err = p.communicate()
        stdout += out.decode()
        stderr += err.decode()

        motifs = []

        if os.path.exists(outfile):
            motifs = self.parse(outfile)
            for m in motifs:
                m.id = "{0}_{1}".format(self.name, m.id)
        else:
            stdout += "\nMotif file {0} not found!\n".format(outfile)
            stderr += "\nMotif file {0} not found!\n".format(outfile)

        return motifs, stdout, stderr

    def parse(self, fname):
        """
        Convert RPMCMC output to motifs

        Parameters
        ----------
        fname : str
            File containing RPMCMC output.

        Returns
        -------
        motifs : list
            List of Motif instances.
        """
        motifs = []
        pfm = []
        name = ""
        for line in open(fname):
            line = line.strip()
            if line.startswith("PFM"):
                continue
            if line.startswith("Motif"):
                if len(pfm) > 0:
                    motif = Motif(pfm)
                    motif.id = name
                    motifs.append(motif)
                name = line
                pfm = []
            else:
                if line != ("A C G T"):
                    row = line.split(" ")
                    if len(row) == 4:
                        row = [float(x) for x in row]
                        pfm.append(row)

        motif = Motif(pfm)
        motif.id = name
        motifs.append(motif)

        return motifs
