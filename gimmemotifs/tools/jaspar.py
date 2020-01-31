from .motifprogram import MotifProgram
import os

from gimmemotifs.motif import read_motifs


class Jaspar(MotifProgram):
    def __init__(self):
        self.name = "JASPAR"
        self.cmd = "/bin/false"
        self.use_width = False

    def _parse_params(self, params=None):
        """
        Parse parameters.

        Combine default and user-defined parameters.
        """
        prm = self.default_params.copy()
        if params is not None:
            prm.update(params)

        # Absolute path, just to be sure
        prm["background"] = os.path.abspath(prm["background"])

        return prm

    def _run_program(self, bin, fastafile, params=None):
        """
        Get enriched JASPAR motifs in a FASTA file.

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
        fname = os.path.join(self.config.get_motif_dir(), "JASPAR2010_vertebrate.pwm")
        motifs = read_motifs(fname, fmt="pwm")

        for motif in motifs:
            motif.id = "JASPAR_%s" % motif.id
        return motifs, "", ""
