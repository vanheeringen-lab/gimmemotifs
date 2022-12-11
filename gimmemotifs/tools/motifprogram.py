import os
from tempfile import mkdtemp

from gimmemotifs.config import MotifConfig
from gimmemotifs.motif import read_motifs


class MotifProgram(object):

    """Motif program base class."""

    config = MotifConfig()
    default_params = {}
    name = None
    tmpdir = None

    def _parse_params(self, params=None, needs_background=False):
        """
        Parse parameters.

        Combine default and user-defined parameters.
        """
        prm = self.default_params.copy()
        if params is not None:
            prm.update(params)

        # Background file is essential!
        if "background" in prm:
            # Absolute path, just to be sure
            prm["background"] = os.path.abspath(prm["background"])
        elif needs_background:
            raise ValueError("Background file needed!")

        return prm

    def _read_and_label_motifs(self, outfile, stdout, stderr, fmt="meme"):
        """Read output motifs and label with program name"""
        if not os.path.exists(outfile):
            stdout += f"\nMotif file {outfile} not found!\n"
            stderr += f"\nMotif file {outfile} not found!\n"
            return [], stdout, stderr

        motifs = read_motifs(outfile, fmt)
        for m in motifs:
            m.id = f"{self.name}_{m.id}"
        return motifs, stdout, stderr

    def bin(self):
        """
        Get the command used to run the tool.

        Returns
        -------
        command : str
            The tool system command.
        """
        return self.config.bin(self.name)

    def dir(self):
        """
        Get the installation directory of the tool.

        Returns
        -------
        dir : str
            The tool directory.
        """
        return self.config.dir(self.name)

    def is_configured(self):
        """
        Check if the tool is configured.

        Returns
        -------
        is_configured : bool
            True if the tool is configured.
        """
        return self.config.is_configured(self.name)

    def is_installed(self):
        """
        Check if the tool is installed.

        Returns
        -------
        is_installed : bool
            True if the tool is installed.
        """
        return self.is_configured() and os.access(self.bin(), os.X_OK)

    def _run_program(self, param, fastafile, params):
        raise NotImplementedError()

    def run(self, fastafile, params=None, tmp=None):
        """
        Run the tool and predict motifs from a FASTA file.

        Parameters
        ----------
        fastafile : str
            Name of the FASTA input file.

        params : dict, optional
            Optional parameters. For some of the tools required parameters
            are passed using this dictionary.

        tmp : str, optional
            Directory to use for creation of temporary files.

        Returns
        -------
        motifs : list of Motif instances
            The predicted motifs.

        stdout : str
            Standard out of the tool.

        stderr : str
            Standard error of the tool.
        """
        if not self.is_configured():
            raise ValueError(f"{self.name} is not configured")

        if not self.is_installed():
            raise ValueError(
                f"{self.name} is not installed or not correctly configured"
            )

        self.tmpdir = mkdtemp(prefix=f"{self.name}.", dir=tmp)
        fastafile = os.path.abspath(fastafile)

        try:
            return self._run_program(self.bin(), fastafile, params)
        except KeyboardInterrupt:
            return [], "Killed", "Killed"
