import platform
import tempfile
import os
import pytest
from gimmemotifs.tools import __tools__, get_tool
from gimmemotifs.comparison import MotifComparer
from gimmemotifs.motif import motif_from_consensus

data_dir = "test/data/motifprogram"
fa = os.path.join("test/data/denovo/input.fa")
bg_fa = os.path.join("test/data/denovo/random.fa")


def ap1_included(motifs):
    ap1 = motif_from_consensus("TGASTCA")
    mc = MotifComparer()
    match = mc.get_closest_match(ap1, motifs, metric="seqcor")
    print(match)
    if match["TGASTCA"][1][0] >= 0.75:
        return True
    return False


@pytest.mark.parametrize("tool_name", __tools__)
def test_tool(tool_name):
    """Test motif prediction tools."""
    params = {"background": bg_fa, "organism": "hg38", "width": 7}
    print(__tools__)
    if tool_name in [
        "yamda",  # not installable vio bioconda
        "rpmcmc",  # not installable via bioconda
        "gadem",  # sometimes crashes on invalid pointer
        "jaspar",
        "xxmotif",  # takes too long
        "trawler",  # unpredictable, sometimes doesn't find the motif
        "weeder",  # doesn't work at the moment
        "posmo",  # motif doesn't predictably look like AP1
        "dreme",  # current dreme in bioconda is broken
    ]:
        return

    if platform.system() == "Darwin":
        # No support for osx
        if tool_name in ["amd", "hms", "improbizer", "motifsampler", "dinamo"]:
            return

    t = get_tool(tool_name)
    print("Testing {}...".format(t))

    (motifs, stderr, stdout) = t.run(fa, params)
    print(motifs)
    print(stderr)
    print(stdout)
    assert ap1_included(motifs)
