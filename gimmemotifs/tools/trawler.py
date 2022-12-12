import glob
import os
import shutil
from subprocess import PIPE, Popen
from tempfile import NamedTemporaryFile

from gimmemotifs.motif import read_motifs

from .motifprogram import MotifProgram


class Trawler(MotifProgram):

    """
    Predict motifs using Trawler.

    Reference: Ettwiller, 2010; PMID: 17589518
    """

    def __init__(self):
        self.name = "Trawler"
        self.cmd = "trawler"
        self.use_width = False
        self.default_params = {"single": False, "background": None}

    def _parse_params(self, params=None):
        """
        Parse parameters.

        Combine default and user-defined parameters.
        """
        prm = super()._parse_params(params)
        prm["strand"] = "double"
        if prm["single"]:
            prm["strand"] = "single"

        return prm

    def _run_program(self, bin, fastafile, params=None):
        """
        Run Trawler and predict motifs from a FASTA file.

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

        tmp = NamedTemporaryFile(mode="w", dir=self.tmpdir, delete=False)
        shutil.copy(fastafile, tmp.name)
        fastafile = tmp.name

        current_path = os.getcwd()
        os.chdir(self.dir())

        motifs = []
        stdout = ""
        stderr = ""
        for wildcard in [0, 1, 2]:
            cmd = (
                f"{bin} -sample {fastafile} -background {params['background']} "
                f"-directory {self.tmpdir} -strand {params['strand']} "
                f"-wildcard {wildcard}"
            )

            p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
            out, err = p.communicate()
            stdout += out.decode()
            stderr += err.decode()

            os.chdir(current_path)
            pfmfiles = glob.glob(f"{self.tmpdir}/tmp*/result/*pwm")
            if len(pfmfiles) > 0:
                out_file = pfmfiles[0]
                stdout += f"\nOutfile: {out_file}"

                my_motifs = []
                if os.path.exists(out_file):
                    my_motifs = read_motifs(out_file, fmt="pwm")
                    for m in motifs:
                        m.id = f"{self.name}_{m.id}"
                    stdout += f"\nTrawler: {len(motifs)} motifs"

                # remove temporary files
                if os.path.exists(tmp.name):
                    os.unlink(tmp.name)

                for motif in my_motifs:
                    motif.id = f"{self.name}_{wildcard}_{motif.id}"

                motifs += my_motifs
            else:
                stderr += "\nNo outfile found"

        return motifs, stdout, stderr
