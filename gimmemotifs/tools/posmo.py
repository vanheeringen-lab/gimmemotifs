from .motifprogram import MotifProgram
import os
import shutil
from subprocess import Popen, PIPE

from gimmemotifs.motif import Motif


class Posmo(MotifProgram):

    """
    Predict motifs using Posmo.

    Reference:
    """

    def __init__(self):
        self.name = "Posmo"
        self.cmd = "posmo"
        self.use_width = True

    def _run_program(self, bin, fastafile, params=None):
        """
        Run Posmo and predict motifs from a FASTA file.

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

        width = params.get("width", 8)
        basename = "posmo_in.fa"

        new_file = os.path.join(self.tmpdir, basename)
        shutil.copy(fastafile, new_file)

        fastafile = new_file
        # pfmfile = fastafile + ".pwm"

        motifs = []
        current_path = os.getcwd()
        os.chdir(self.tmpdir)
        for n_ones in range(4, min(width, 11), 2):
            x = "1" * n_ones
            outfile = "%s.%s.out" % (fastafile, x)
            cmd = "%s 5000 %s %s 1.6 2.5 %s 200" % (bin, x, fastafile, width)
            p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
            stdout, stderr = p.communicate()
            stdout = stdout.decode()
            stderr = stderr.decode()

            context_file = fastafile.replace(
                basename, "context.%s.%s.txt" % (basename, x)
            )
            cmd = "%s %s %s simi.txt 0.88 10 2 10" % (
                bin.replace("posmo", "clusterwd"),
                context_file,
                outfile,
            )
            p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
            out, err = p.communicate()
            stdout += out.decode()
            stderr += err.decode()

            if os.path.exists(outfile):
                with open(outfile) as f:
                    motifs += self.parse(f, width, n_ones)

        os.chdir(current_path)

        return motifs, stdout, stderr

    def parse(self, fo, width, seed=None):
        """
        Convert Posmo output to motifs

        Parameters
        ----------
        fo : file-like
            File object containing Posmo output.

        Returns
        -------
        motifs : list
            List of Motif instances.
        """
        motifs = []

        lines = [fo.readline() for x in range(6)]
        while lines[0]:
            matrix = [
                [float(x) for x in line.strip().split("\t")] for line in lines[2:]
            ]
            matrix = [[matrix[x][y] for x in range(4)] for y in range(len(matrix[0]))]
            m = Motif(matrix)
            m.trim(0.1)
            m.id = lines[0].strip().split(" ")[-1]
            motifs.append(m)
            lines = [fo.readline() for x in range(6)]

        for i, motif in enumerate(motifs):
            if seed:
                motif.id = "%s_w%s.%s_%s" % (self.name, width, seed, i + 1)
            else:
                motif.id = "%s_w%s_%s" % (self.name, width, i + 1)
            motif.trim(0.25)

        return motifs
