import os
from tempfile import TemporaryDirectory

import pytest

from gimmemotifs.cli import cli
from gimmemotifs.motif import read_motifs

travis = "TRAVIS" in os.environ and os.environ["TRAVIS"] == "true"


# @pytest.mark.skipif(travis, reason="Skip CPU-intensive tests")
@pytest.mark.parametrize(
    "motif_argument",
    [
        pytest.param(["--known"], id="known"),
        pytest.param(["--denovo"], id="denovo"),
        pytest.param([], id="default"),
    ],
)
def test_gimme_motifs(motif_argument):
    with TemporaryDirectory() as d:
        print(d)
        cli(
            ["motifs", "test/data/denovo/input.fa", d] +  # test/data/cli/Gm12878.CTCF.top200.fa
            ["-p", "test/data/cli/motifs.pfm"] +
            ["-g", "test/data/background/genome.fa"] +
            ["-a", "small", "-t", "MEME", "--nogc", "-N", "1"] +
            motif_argument
        )

        assert 1 == 1


# @pytest.mark.skipif(travis, reason="Skip CPU-intensive tests")
def test_gimme_maelstrom():
    with TemporaryDirectory() as d:
        print(d)
        cli(
            ["maelstrom", "test/data/maelstrom/input_table.txt"] +  # test/data/maelstrom/input.table.txt
            ["test/data/background/genome.fa", d] +
            ["--nogc", "-m", "RF"]  # "--no-filter",
        )

        assert 1 == 1


@pytest.mark.parametrize(
    "arguments",
    [
        ["-c", "0.8"],
        ["-t", "-g", "test/data/genomes/hg38sample.fa"],
        ["-T"],
        ["-b", "-g" "test/data/genomes/hg38sample.fa"],
        ["-z", "--gc", "-g", "test/data/genomes/hg38sample.fa"],
    ],
)
def test_gimme_scan(arguments):
    cli(
        ["scan", "test/data/pwmscan/10promoters.fa", "-p", "test/data/pwmscan/TATA.pwm"]
        + arguments
    )

    assert 1 == 1


def test_gimme_logo():
    motif_name = "MA0103.3_ZEB1"
    cli(["logo", "-p", "test/data/cli/motifs.pfm", "-i", motif_name])

    assert os.path.exists(f"{motif_name}.png")
    if os.path.exists(f"{motif_name}.png"):
        os.unlink(f"{motif_name}.png")


def test_gimme_match(capsys):
    pfm = "test/data/cli/test_motif.pfm"
    for db_args, out in [
        ([], "GM.5.0.p53.0001"),
        (["-d", "JASPAR2020_vertebrates"], "MA0106.3_TP53"),
    ]:

        cli(["match", pfm] + db_args)
        captured = capsys.readouterr()
        match = captured.out.strip().split("\n")[-1].split("\t")[1]
        assert match == out


def test_gimme_cluster():
    pfm = "test/data/cli/cluster.pfm"

    with TemporaryDirectory() as d:
        cli(["cluster", pfm, d, "-t", "0.99"])
        out_pfm = os.path.join(d, "clustered_motifs.pfm")
        assert os.path.exists(out_pfm)
        cons = ["ACCGTTAACsGy", "ATGACkyA", "TTGCGnAA"]
        assert cons == sorted([m.to_consensus() for m in read_motifs(out_pfm)])
