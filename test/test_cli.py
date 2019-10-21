import os
from tempfile import TemporaryDirectory

import pytest

from gimmemotifs.cli import cli

travis = "TRAVIS" in os.environ and os.environ["TRAVIS"] == "true"


@pytest.mark.skipif(travis, reason="Skip CPU-intensive tests")
@pytest.mark.parametrize(
    "denovo_known",
    [
        ["-p", "test/data/cli/motifs.pfm", "-t", "MDmodule"],
        ["--denovo", "-t", "MDmodule"],
        ["--known", "-p", "test/data/cli/motifs.pfm"],
    ],
)
def test_gimme_motifs(denovo_known):
    with TemporaryDirectory() as d:
        print(d)
        cli(
            ["motifs", "test/data/cli/Gm12878.CTCF.top200.fa", d, "-g", "hg19"]
            + denovo_known
        )

        assert 1 == 1


@pytest.mark.skipif(travis, reason="Skip CPU-intensive tests")
def test_gimme_maelstrom():
    with TemporaryDirectory() as d:
        print(d)
        cli(["maelstrom", "test/data/maelstrom/input.table.txt", "mm10", d])

        assert 1 == 1


@pytest.mark.skipif(travis, reason="Skip CPU-intensive tests")
@pytest.mark.parametrize(
    "arguments",
    [
        ["-c", "0.8"],
        ["-t", "-g", "hg38"],
        ["-T"],
        ["-b", "-g" "hg38"],
        ["-z", "--gc", "-g", "hg38"],
    ],
)
def test_gimme_scan(arguments):
    cli(
        ["scan", "test/data/pwmscan/promoters.fa", "-p", "test/data/pwmscan/TATA.pwm"]
        + arguments
    )

    assert 1 == 1

def test_gimme_logo():
    motif_name = "MA0103.3_ZEB1"
    cli(
        ["logo", "test/data/cli/motifs.pfm", "-i", motif_name]
    )

    assert os.path.exists(f"{motif_name}.png")
    if os.path.exists(f"{motif_name}.png"):
        os.unlink(f"{motif_name}.png")
