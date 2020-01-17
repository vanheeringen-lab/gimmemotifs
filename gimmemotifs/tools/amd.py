from .motifprogram import MotifProgram
import os
import re
import shutil
from subprocess import Popen, PIPE

from gimmemotifs.motif import Motif


class Amd(MotifProgram):

    """
    Predict motifs using AMD.

    Reference:
    """

    def __init__(self):
        self.name = "AMD"
        self.cmd = "AMD.bin"
        self.use_width = False
        self.default_params = {"background": None}

    def _run_program(self, bin, fastafile, params=None):
        """
        Run AMD and predict motifs from a FASTA file.

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
        params = self._parse_params(params, needs_background=True)

        fgfile = os.path.join(self.tmpdir, "AMD.in.fa")
        outfile = fgfile + ".Matrix"
        shutil.copy(fastafile, fgfile)

        current_path = os.getcwd()
        os.chdir(self.tmpdir)

        stdout = ""
        stderr = ""

        cmd = "%s -F %s -B %s" % (bin, fgfile, params["background"])
        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()
        stdout += out.decode()
        stderr += err.decode()

        os.chdir(current_path)
        motifs = []
        if os.path.exists(outfile):
            f = open(outfile)
            motifs = self.parse(f)
            f.close()

        return motifs, stdout, stderr

    def parse(self, fo):
        """
        Convert AMD output to motifs

        Parameters
        ----------
        fo : file-like
            File object containing AMD output.

        Returns
        -------
        motifs : list
            List of Motif instances.
        """
        motifs = []

        # 160:  112  CACGTGC      7.25   chr14:32308489-32308689
        p = re.compile(r"\d+\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)")
        wm = []
        name = ""
        for line in fo.readlines():
            if line.startswith("Motif") and line.strip().endswith(":"):
                if name:
                    motifs.append(Motif(wm))
                    motifs[-1].id = name
                    name = ""
                    wm = []
                name = "%s_%s" % (self.name, line.split(":")[0])
            else:
                m = p.search(line)
                if m:
                    wm.append([float(m.group(x)) for x in range(1, 5)])
        motifs.append(Motif(wm))
        motifs[-1].id = name

        return motifs
