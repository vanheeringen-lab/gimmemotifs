from .motifprogram import MotifProgram
import os
from subprocess import Popen, PIPE

from gimmemotifs.motif import read_motifs


class Dreme(MotifProgram):

    """
    Predict motifs using DREME.

    Reference: Bailey, 2011, https://doi.org/10.1093/bioinformatics/btr261
    """

    def __init__(self):
        self.name = "DREME"
        self.cmd = "dreme-py3"
        self.use_width = True

    def _run_program(self, bin, fastafile, params=None):
        """
        DREME Run  and predict motifs from a FASTA file.

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
        default_params = {"single": False, "number": 10}
        if params is not None:
            default_params.update(params)

        outfile = os.path.join(self.tmpdir, "dreme.txt")

        strand = " -norc "
        number = default_params["number"]

        cmd = [bin, "-p", fastafile, "-m", "%s" % number, "-oc", self.tmpdir]
        if default_params["background"]:
            cmd += ["-n", default_params["background"]]
        if default_params["single"]:
            cmd.append(strand)

        p = Popen(cmd, bufsize=1, stderr=PIPE, stdout=PIPE)
        stdout, stderr = p.communicate()

        motifs = read_motifs(outfile, fmt="meme")
        for motif in motifs:
            motif.id = self.name + "_" + motif.id

        return motifs, stdout, stderr
