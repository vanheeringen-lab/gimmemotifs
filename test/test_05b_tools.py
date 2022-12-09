import platform
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
    blacklist = {
        "dreme": "not installable via bioconda/source",
        "gadem": "sometimes crashes on invalid pointer",
        "jaspar": "# TODO: not configured",  # TODO
        "posmo": "motif doesn't predictably look like AP1",
        "rpmcmc": "not installable via bioconda/source",
        "trawler": "unpredictable, sometimes doesn't find the motif",
        "xxmotif": "takes too long",
        "yamda": "# TODO: not configured",  # TODO
    }
    if tool_name in blacklist:
        pytest.skip(blacklist[tool_name])

    if platform.system() == "Darwin":
        if tool_name in ["amd", "dinamo", "hms", "improbizer", "motifsampler"]:
            pytest.skip("No supported for osx")

    print("Tool class:", __tools__[tool_name])
    t = get_tool(tool_name)
    print(f"Testing {t}...")

    params = {"background": bg_fa, "organism": "hg38", "width": 7}
    (motifs, stderr, stdout) = t.run(fa, params)
    print(motifs)
    print(stderr)
    print(stdout)
    assert ap1_included(motifs)
