from .motifprogram import MotifProgram
import os
import shutil
from subprocess import Popen, PIPE

from gimmemotifs.motif import Motif


class Gadem(MotifProgram):

    """
    Predict motifs using GADEM.

    Reference:
    """

    def __init__(self):
        self.name = "GADEM"
        self.cmd = "gadem"
        self.use_width = False

    def _run_program(self, bin, fastafile, params=None):
        """
        Run GADEM and predict motifs from a FASTA file.

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
        default_params = {}
        if params is not None:
            default_params.update(params)

        new_file = os.path.join(self.tmpdir, "gadem_in.fa")
        shutil.copy(fastafile, new_file)

        fastafile = new_file
        pfmfile = fastafile + ".pwm"
        outfile = fastafile + ".out"

        current_path = os.getcwd()
        os.chdir(self.tmpdir)
        cmd = "%s -fseq %s -fpwm %s -fout %s" % (bin, fastafile, pfmfile, outfile)
        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        stdout, stderr = p.communicate()

        motifs = []
        if os.path.exists(pfmfile):
            with open(pfmfile) as f:
                motifs = self.parse(f)

        os.chdir(current_path)

        return motifs, stdout, stderr

    def parse(self, fo):
        """
        Convert GADEM output to motifs

        Parameters
        ----------
        fo : file-like
            File object containing GADEM output.

        Returns
        -------
        motifs : list
            List of Motif instances.
        """
        motifs = []
        nucs = {"A": 0, "C": 1, "G": 2, "T": 3}

        lines = fo.readlines()
        for i in range(0, len(lines), 5):
            align = []
            pwm = []
            pfm = []
            m_id = ""
            line = lines[i].strip()
            m_id = line[1:]
            number = m_id.split("_")[0][1:]
            if os.path.exists("%s.seq" % number):
                with open("%s.seq" % number) as f:
                    for line in f:
                        if "x" not in line and "n" not in line:
                            line = line.strip().upper()
                            align.append(line)
                            if not pfm:
                                pfm = [[0 for x in range(4)] for x in range(len(line))]
                            for p in range(len(line)):
                                pfm[p][nucs[line[p]]] += 1

            m = [
                line.strip().split(" ")[1].split("\t") for line in lines[i + 1 : i + 5]
            ]

            pwm = [[float(m[x][y]) for x in range(4)] for y in range(len(m[0]))]

            motifs.append(Motif(pwm))
            motifs[-1].id = "{}_{}".format(self.name, m_id)
            # motifs[-1].pwm = pwm
            if align:
                motifs[-1].pfm = pfm
                motifs[-1].align = align

        return motifs
