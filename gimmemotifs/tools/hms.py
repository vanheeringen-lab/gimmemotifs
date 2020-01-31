from .motifprogram import MotifProgram
import os
import shutil
from subprocess import Popen, PIPE

from gimmemotifs.fasta import Fasta
from gimmemotifs.motif import Motif


class Hms(MotifProgram):

    """
    Predict motifs using HMS.

    Reference:
    """

    def __init__(self):
        self.name = "HMS"
        self.cmd = "hms"
        self.use_width = True
        self.default_params = {"background": None}

    def _prepare_files(self, fastafile):

        hmsdir = self.dir()
        thetas = ["theta%s.txt" % i for i in [0, 1, 2, 3]]
        for t in thetas:
            shutil.copy(os.path.join(hmsdir, t), self.tmpdir)

        summitfile = os.path.join(self.tmpdir, "HMS.in.summits.txt")
        outfile = os.path.join(self.tmpdir, "thetafinal.txt")
        fgfile = os.path.join(self.tmpdir, "HMS.in.fa")

        shutil.copy(fastafile, fgfile)
        fa = Fasta(fgfile)
        with open(summitfile, "w") as out:
            for seq in fa.seqs:
                out.write("%s\n" % (len(seq) / 2))
        return fgfile, summitfile, outfile

    def _run_program(self, bin, fastafile, params=None):
        """
        Run HMS and predict motifs from a FASTA file.

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

        default_params = {"width": 10}
        if params is not None:
            default_params.update(params)

        fgfile, summitfile, outfile = self._prepare_files(fastafile)

        current_path = os.getcwd()
        os.chdir(self.tmpdir)

        cmd = (
            "{} -i {} -w {} -dna 4 -iteration 50 -chain 20 -seqprop -0.1 "
            "-strand 2 -peaklocation {} -t_dof 3 -dep 2"
        ).format(bin, fgfile, params["width"], summitfile)

        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        stdout, stderr = p.communicate()

        os.chdir(current_path)

        motifs = []
        if os.path.exists(outfile):
            with open(outfile) as f:
                motifs = self.parse(f)
                for i, m in enumerate(motifs):
                    m.id = "HMS_w{}_{}".format(params["width"], i + 1)

        return motifs, stdout, stderr

    def parse(self, fo):
        """
        Convert HMS output to motifs

        Parameters
        ----------
        fo : file-like
            File object containing HMS output.

        Returns
        -------
        motifs : list
            List of Motif instances.
        """
        motifs = []
        m = [[float(x) for x in fo.readline().strip().split(" ")] for i in range(4)]
        matrix = [[m[0][i], m[1][i], m[2][i], m[3][i]] for i in range(len(m[0]))]
        motifs = [Motif(matrix)]
        motifs[-1].id = self.name

        return motifs
