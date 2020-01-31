from .motifprogram import MotifProgram
import os
import re
import shutil
from subprocess import Popen, PIPE

from gimmemotifs.motif import Motif


class MDmodule(MotifProgram):

    """
    Predict motifs using MDmodule.

    Reference:
    """

    def __init__(self):
        self.name = "MDmodule"
        self.cmd = "MDmodule"
        self.use_width = True

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
        Run MDmodule and predict motifs from a FASTA file.

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
        default_params = {"width": 10, "number": 10}
        if params is not None:
            default_params.update(params)

        new_file = os.path.join(self.tmpdir, "mdmodule_in.fa")
        shutil.copy(fastafile, new_file)

        fastafile = new_file
        pfmfile = fastafile + ".out"

        width = default_params["width"]
        number = default_params["number"]

        current_path = os.getcwd()
        os.chdir(self.tmpdir)
        cmd = "%s -i %s -a 1 -o %s -w %s -t 100 -r %s" % (
            bin,
            fastafile,
            pfmfile,
            width,
            number,
        )
        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        stdout, stderr = p.communicate()

        stdout = "cmd: {}\n".format(cmd) + stdout.decode()

        motifs = []
        if os.path.exists(pfmfile):
            with open(pfmfile) as f:
                motifs = self.parse(f)

        os.chdir(current_path)

        for motif in motifs:
            motif.id = "%s_%s" % (self.name, motif.id)

        return motifs, stdout, stderr

    def parse(self, fo):
        """
        Convert MDmodule output to motifs

        Parameters
        ----------
        fo : file-like
            File object containing MDmodule output.

        Returns
        -------
        motifs : list
            List of Motif instances.
        """
        motifs = []
        nucs = {"A": 0, "C": 1, "G": 2, "T": 3}
        p = re.compile(r"(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)")
        pf = re.compile(r">.+\s+[bf]\d+\s+(\w+)")

        pwm = []
        pfm = []
        align = []
        m_id = ""
        for line in fo.readlines():
            if line.startswith("Motif"):
                if m_id:
                    motifs.append(Motif())
                    motifs[-1].id = m_id
                    motifs[-1].pwm = pwm
                    motifs[-1].pfm = pfm
                    motifs[-1].align = align
                    pwm = []
                    pfm = []
                    align = []
                m_id = line.split("\t")[0]
            else:
                m = p.search(line)
                if m:
                    pwm.append([float(m.group(x)) / 100 for x in [2, 3, 4, 5]])
                m = pf.search(line)
                if m:
                    if not pfm:
                        pfm = [[0 for x in range(4)] for x in range(len(m.group(1)))]
                    for i in range(len(m.group(1))):
                        pfm[i][nucs[m.group(1)[i]]] += 1

                    align.append(m.group(1))

        if pwm:
            motifs.append(Motif())
            motifs[-1].id = m_id
            motifs[-1].pwm = pwm
            motifs[-1].pfm = pfm
            motifs[-1].align = align

        return motifs
