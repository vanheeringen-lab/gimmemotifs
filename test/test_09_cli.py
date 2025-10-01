import os
import sys
from tempfile import TemporaryDirectory
from unittest.mock import patch

import pytest

from gimmemotifs.cli import gimme
from gimmemotifs.motif import read_motifs

travis = "TRAVIS" in os.environ and os.environ["TRAVIS"] == "true"


# @pytest.mark.skipif(travis, reason="Skip CPU-intensive tests")
@pytest.mark.parametrize(
    "args",
    [
        pytest.param(["--known"], id="known"),
        pytest.param(["--denovo"], id="denovo"),
        pytest.param([], id="default"),
    ],
)
def test_gimme_motifs(args):
    input = "test/data/denovo/input.fa"
    with (TemporaryDirectory() as outdir):
        cmd = ["gimme", "motifs", input, outdir, *args]
        cmd += ["-p", "test/data/cli/motifs.pfm"]
        cmd += ["-g", "test/data/background/genome.fa"]
        cmd += ["-a", "small", "-t", "MEME", "--nogc", "-N", "1"]
        with patch.object(sys, 'argv', cmd):
            print(sys.argv)
            gimme()

        observed = sorted(os.listdir(outdir))
    if "--known" in args:
        expected = [
            'generated_background.gc.fa',
            'gimme.motifs.html',
            'gimme.motifs.redundant.html',
            'gimme.roc.report.txt',
            'logos',
            'motif_scan_results',
        ]
    else:
        expected = [
            'combined.motif2factors.txt',
            'combined.pfm',
            'generated_background.gc.fa',
            'gimme.clustereds.html',
            'gimme.denovo.html',
            'gimme.denovo.pfm',
            'gimme.motifs.html',
            'gimme.motifs.redundant.html',
            'gimme.roc.report.txt',
            'gimmemotifs.log',
            'images',
            'logos',
            'motif_scan_results',
            'params.txt',
            'stats.gc.txt',
        ]
    assert observed == expected


@pytest.mark.skipif(travis, reason="Can cause a seg fault on Travis")
def test_gimme_maelstrom():
    input = "test/data/maelstrom/input_table.txt"
    genome = "test/data/background/genome.fa"
    with TemporaryDirectory() as outdir:
        cmd = [
                "gimme",
                "maelstrom",
                input,
                genome,
                outdir,
                "--nogc",
                "-m RF",
                "-s 123",
            ]
        with patch.object(sys, 'argv', cmd):
            print(sys.argv)
            gimme()
        observed = sorted(os.listdir(outdir))
    expected = [
        'activity.rf.score.out.txt',
        'gimme.vertebrate.v5.0.motif2factors.txt',
        'gimme.vertebrate.v5.0.pfm',
        'input.table.txt',
        'motif.count.txt.gz',
        'motif.freq.txt',
        'motif.nr.count.txt.gz',
        'motif.nr.score.txt.gz',
        'motif.score.txt.gz',
        'nonredundant.motifs.motif2factors.txt',
        'nonredundant.motifs.pfm',
    ]
    assert observed == expected


@pytest.mark.parametrize(
    "args",
    [
        ["-c", "0.8"],
        ["-t", "-g", "test/data/genomes/hg38sample.fa"],
        ["-T"],
        ["-b", "-g" "test/data/genomes/hg38sample.fa"],
        ["-z", "--gc", "-g", "test/data/genomes/hg38sample.fa"],
        ["-s", "123", "-g", "test/data/genomes/hg38sample.fa"],
    ],
)
def test_gimme_scan(args):
    input = "test/data/pwmscan/10promoters.fa"
    cmd = ["gimme", "scan", input, "-p", "test/data/pwmscan/TATA.pwm", *args]
    with patch.object(sys, 'argv', cmd):
        gimme()

    assert 1 == 1


def test_gimme_logo():
    motif_name = "MA0103.3_ZEB1"
    cmd = ["gimme", "logo", "-p", "test/data/cli/motifs.pfm", "-i", motif_name]
    with patch.object(sys, 'argv', cmd):
        gimme()

    assert os.path.exists(f"{motif_name}.png")
    if os.path.exists(f"{motif_name}.png"):
        os.unlink(f"{motif_name}.png")


@pytest.mark.parametrize(
    ("db_args", "out"),
    [
        ([], "GM.5.0.p53.0001"),
        (["-d", "JASPAR2020_vertebrates"], "MA0106.3_TP53"),
    ],
)
def test_gimme_match(db_args, out, capsys):
    cmd = ["gimme", "match", "test/data/cli/test_motif.pfm", *db_args]
    with patch.object(sys, 'argv', cmd):
        gimme()
        captured = capsys.readouterr()
    match = captured.out.strip().split("\n")[-1].split("\t")[1]
    assert match == out


def test_gimme_cluster():
    with TemporaryDirectory() as outdir:
        cmd = ["gimme", "cluster", "test/data/cli/cluster.pfm", outdir, "-t", "0.99"]
        with patch.object(sys, 'argv', cmd):
            gimme()
        out_pfm = os.path.join(outdir, "clustered_motifs.pfm")
        assert os.path.exists(out_pfm)
        cons = ["ACCGTTAACsGy", "ATGACkyA", "TTGCGnAA"]
        assert cons == sorted([m.to_consensus() for m in read_motifs(out_pfm)])
