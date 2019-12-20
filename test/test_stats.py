import os

import pytest

from gimmemotifs.stats import calc_stats
from gimmemotifs.motif import read_motifs
from gimmemotifs import rocmetrics


data_dir = "test/data/stats"
genome = "test/data/background/genome.fa"
motifs = os.path.join(data_dir, "motifs.pwm")
fg_fa = os.path.join(data_dir, "p73.fa")
bg_fa = os.path.join(data_dir, "random.w200.fa")
fg_table = os.path.join(data_dir, "p73.scores.txt")
bg_table = os.path.join(data_dir, "random.scores.txt")


@pytest.fixture
def stat_functions():
    return [
        "recall_at_fdr",
        "fraction_fpr",
        "score_at_fpr",
        "enr_at_fpr",
        "max_enrichment",
        "mncp",
        "roc_auc",
        "roc_auc_xlim",
        "max_fmeasure",
        "ks_pvalue",
        "ks_significance",
    ]


@pytest.mark.parametrize(
    "kwargs",
    [
        {
            "motifs": motifs,
            "fg_file": fg_fa,
            "bg_file": bg_fa,
            "genome": genome,
            "gc": True,
            "zscore": True,
        },
        {
            "motifs": motifs,
            "fg_file": fg_fa,
            "bg_file": bg_fa,
            "genome": genome,
            "gc": False,
            "zscore": False,
        },
        {"motifs": motifs, "fg_table": fg_table, "bg_table": bg_table,},
    ],
)
def test1_stats(kwargs, stat_functions):
    """ Calculate motif statistics """
    for ncpus in [1, 2]:
        kwargs["ncpus"] = ncpus
        stats = calc_stats(**kwargs)

        for f in stat_functions:
            if "fg_table" not in kwargs or getattr(rocmetrics, f).input_type != "pos":
                print(f, fg_table, getattr(rocmetrics, f).input_type)
                assert f in list(stats.values())[0]

        # Two motifs
        assert 2 == len(stats)

        m1 = "T-box_M1713_1.01_CTAGGTGTGAA"  # not enriched
        m2 = "p53_Average_8_CATGyCnGGrCATGy"  # highly enriched

        assert stats[m1]["roc_auc"] < 0.9
        assert stats[m2]["roc_auc"] > 0.5

        assert stats[m1]["recall_at_fdr"] == 0.0
        assert stats[m2]["recall_at_fdr"] > 0.8

        if "fg_table" not in kwargs:
            assert stats[m1]["ks_pvalue"] > 0.01
            assert stats[m2]["ks_pvalue"] < 0.001

        assert stats[m1]["phyper_at_fpr"] > 0.1
        assert stats[m2]["phyper_at_fpr"] < 1e-13


@pytest.mark.parametrize(
    "kwargs",
    [
        {
            "motifs": motifs,
            "fg_file": fg_fa,
            "bg_file": bg_fa,
            "stats": ["roc_auc"],
            "gc": False,
            "zscore": False,
        },
    ],
)
def test_one_statistic(kwargs):
    # Only calculate specific statistic
    stats = calc_stats(**kwargs)

    assert 1 == len(list(stats.values())[0])

    m1 = "T-box_M1713_1.01_CTAGGTGTGAA"  # not enriched
    m2 = "p53_Average_8_CATGyCnGGrCATGy"  # highly enriched
    assert stats[m1]["roc_auc"] < 0.9
    assert stats[m2]["roc_auc"] > 0.5


@pytest.mark.parametrize(
    "kwargs",
    [
        {
            "motifs": motifs,
            "fg_file": fg_fa,
            "bg_file": bg_fa,
            "stats": ["roc_auc"],
            "gc": False,
            "zscore": False,
        },
    ],
)
def test2_stats_single_motif(kwargs):
    """ Calculate motif statistics """

    m_id = "p53_Average_8_CATGyCnGGrCATGy"

    motifs = read_motifs(kwargs["motifs"])
    motif = [m for m in motifs if str(m) == m_id][0]
    kwargs["motifs"] = motif

    stats = calc_stats(**kwargs)
    assert stats[m_id]["roc_auc"] > 0.9
