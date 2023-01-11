import os
from subprocess import PIPE, Popen
from tempfile import NamedTemporaryFile

from gimmemotifs.motif import Motif

from .motifprogram import MotifProgram


class MotifSampler(MotifProgram):

    """
    Predict motifs using MotifSampler.

    Reference:
    """

    def __init__(self):
        self.name = "MotifSampler"
        self.cmd = "MotifSampler"
        self.use_width = True
        self.default_params = {
            "width": 10,
            "background_model": "",
            "single": False,
            "number": 10,
        }

    def _parse_params(self, params=None):
        """
        Parse parameters.

        Combine default and user-defined parameters.
        """
        prm = super()._parse_params(params)

        if prm["background_model"]:
            # Absolute path, just to be sure
            prm["background_model"] = os.path.abspath(prm["background_model"])
        else:
            if prm.get("organism", None):
                prm["background_model"] = os.path.join(
                    self.config.get_bg_dir(),
                    f"{prm['organism']}.MotifSampler.bg",
                )
            else:
                raise Exception(f"No background specified for {self.name}")

        prm["strand"] = 1
        if prm["single"]:
            prm["strand"] = 0

        tmp = NamedTemporaryFile(dir=self.tmpdir)
        prm["pfmfile"] = tmp.name

        tmp2 = NamedTemporaryFile(dir=self.tmpdir)
        prm["outfile"] = tmp2.name

        return prm

    def _run_program(self, bin, fastafile, params=None):
        """
        Run MotifSampler and predict motifs from a FASTA file.

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
        # TODO: test organism

        cmd = (
            f"{bin} -f {fastafile} -b {params['background_model']} "
            f"-m {params['pfmfile']} -w {params['width']} -n {params['number']} "
            f"-o {params['outfile']} -s {params['strand']}"
        )

        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        stdout, stderr = p.communicate()

        motifs = []
        if os.path.exists(params["outfile"]):
            with open(params["outfile"]) as f:
                motifs = self.parse_out(f)

        for motif in motifs:
            motif.id = f"{self.name}_{motif.id}"

        return motifs, stdout, stderr

    def parse(self, fo):
        """
        Convert MotifSampler output to motifs

        Parameters
        ----------
        fo : file-like
            File object containing MotifSampler output.

        Returns
        -------
        motifs : list
            List of Motif instances.
        """
        motifs = []

        ppm = []
        info = {}
        for line in fo.readlines():
            if line.startswith("#"):
                vals = line.strip()[1:].split(" = ")
                if len(vals) > 1:
                    info[vals[0]] = vals[1]
            elif len(line) > 1:
                ppm.append([float(x) for x in line.strip().split("\t")])
            else:
                motifs.append(Motif(ppm=ppm[:]))
                motifs[-1].consensus = info["Consensus"]
                motifs[-1].width = info["W"]
                motifs[-1].id = info["ID"]
                ppm = []

        return motifs

    def parse_out(self, fo):
        """
        Convert MotifSampler output to motifs

        Parameters
        ----------
        fo : file-like
            File object containing MotifSampler output.

        Returns
        -------
        motifs : list
            List of Motif instances.
        """
        motifs = []
        nucs = {"A": 0, "C": 1, "G": 2, "T": 3}
        pseudo = 0.0  # Should be 1/sqrt(# of seqs)
        aligns = {}
        for line in fo.readlines():
            if line.startswith("#"):
                pass
            elif len(line) > 1:
                vals = line.strip().split("\t")
                m_id, site = [
                    x.strip().split(" ")[1].replace('"', "")
                    for x in vals[8].split(";")
                    if x
                ]
                # if vals[6] == "+":
                if site.upper().find("N") == -1:
                    aligns.setdefault(m_id, []).append(site)
                # else:
                #    print site, rc(site)
                #    aligns.setdefault(id, []).append(rc(site))

        for m_id, align in aligns.items():
            # print id, len(align)

            width = len(align[0])
            pfm = [[0 for x in range(4)] for x in range(width)]
            for row in align:
                for i in range(len(row)):
                    pfm[i][nucs[row[i]]] += 1
            total = float(len(align))
            ppm = [[(x + pseudo / 4) / total + (pseudo) for x in row] for row in pfm]
            m = Motif(pfm=pfm, ppm=ppm)
            m.align = align[:]
            m.id = m_id
            motifs.append(m)
        return motifs
