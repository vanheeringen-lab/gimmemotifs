import os
from subprocess import PIPE, Popen

from gimmemotifs.fasta import Fasta
from gimmemotifs.motif import Motif

from .motifprogram import MotifProgram


class ChIPMunk(MotifProgram):

    """
    Predict motifs using ChIPMunk.

    Reference:
    """

    def __init__(self):
        self.name = "ChIPMunk"
        self.cmd = "ChIPMunk.sh"
        self.use_width = True
        self.default_params = {}

    def _run_program(self, bin, fastafile, params=None):
        """
        Run ChIPMunk and predict motifs from a FASTA file.

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
        basename = "munk_in.fa"

        new_file = os.path.join(self.tmpdir, basename)
        out = open(new_file, "w")
        f = Fasta(fastafile)
        for seq in f.seqs:
            header = len(seq) // 2
            out.write(f">{header}\n")
            out.write(f"{seq}\n")
        out.close()

        fastafile = new_file
        outfile = fastafile + ".out"

        current_path = os.getcwd()
        os.chdir(self.dir())

        motifs = []
        # Max recommended by ChIPMunk userguide
        ncpus = 4
        stdout = ""
        stderr = ""
        for zoops_factor in ["oops", 0.0, 0.5, 1.0]:
            cmd = (
                f"{bin} {params.get('width', 8)} {params.get('width', 20)} "
                f"y {zoops_factor} m:{fastafile} 100 10 1 {ncpus} 1>{outfile}"
            )

            p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
            std = p.communicate()
            stdout = stdout + std[0].decode()
            stderr = stderr + std[1].decode()

            if "RuntimeException" in stderr:
                return [], stdout, stderr

            if os.path.exists(outfile):
                with open(outfile) as f:
                    motifs += self.parse(f)

        os.chdir(current_path)

        return motifs, stdout, stderr

    def parse(self, fo):
        """
        Convert ChIPMunk output to motifs

        Parameters
        ----------
        fo : file-like
            File object containing ChIPMunk output.

        Returns
        -------
        motifs : list
            List of Motif instances.
        """
        line = fo.readline()
        if not line:
            return []

        while not line.startswith("A|"):
            line = fo.readline()
        matrix = []
        for _ in range(4):
            matrix.append([float(x) for x in line.strip().split("|")[1].split(" ")])
            line = fo.readline()
        # print matrix
        matrix = [[matrix[x][y] for x in range(4)] for y in range(len(matrix[0]))]
        # print matrix
        m = Motif(matrix)
        m.id = f"ChIPMunk_w{len(m)}"
        return [m]
