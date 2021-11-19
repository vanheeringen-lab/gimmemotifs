from io import StringIO
import os
import pytest
import numpy as np

from gimmemotifs.motif import Motif, read_motifs, motif_from_consensus


@pytest.fixture
def pfm():
    return [
        [3, 3, 3, 3],
        [12, 0, 0, 0],
        [0, 12, 0, 0],
        [0, 0, 12, 0],
        [0, 0, 0, 12],
        [6, 6, 0, 0],
        [0, 6, 6, 0],
        [0, 0, 6, 6],
        [6, 0, 6, 0],
        [6, 0, 0, 6],
        [0, 6, 0, 6],
        [3, 3, 3, 3],
    ]


@pytest.fixture
def my_motif(pfm):
    return Motif(pfm)


@pytest.fixture
def data_dir():
    return "test/data/"


@pytest.fixture
def pfmfile(data_dir):
    return os.path.join(data_dir, "motif/test.pwm")


@pytest.fixture
def pfmfile2(data_dir):
    return os.path.join(data_dir, "pwms/motifs.pwm")


@pytest.fixture
def jaspar(data_dir):
    return os.path.join(data_dir, "pwms/test.jaspar")


def test_motif_properties(pfm):
    m = Motif(pfm)
    assert m.pfm.shape == (12, 4)
    assert m.ppm.shape == (12, 4)
    np.testing.assert_allclose(m.ppm[0], [0.25, 0.25, 0.25, 0.25], atol=0.001)
    np.testing.assert_allclose(m.ppm[1], [1, 0, 0, 0], atol=0.001)
    np.testing.assert_allclose(m.ppm[-2], [0, 0.5, 0, 0.5], atol=0.001)
    assert m.logodds.shape == (12, 4)


def test_motif_length(my_motif):
    assert 12 == len(my_motif)


def test_motif_consensus(my_motif):
    assert 12 == len(my_motif.consensus)
    assert "nACGTmskrwyn" == my_motif.consensus


def test_scores(my_motif):
    np.testing.assert_almost_equal(my_motif.min_score, -45.6396, decimal=4)
    np.testing.assert_almost_equal(my_motif.max_score, 9.7615, decimal=4)


def test_motif_slice(my_motif):
    m_new = my_motif[2:5]

    assert 3 == len(m_new)
    assert "CGT" == m_new.consensus
    np.testing.assert_allclose(m_new.pfm, [[0, 12, 0, 0], [0, 0, 12, 0], [0, 0, 0, 12]])
    np.testing.assert_allclose(
        m_new.ppm, [[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], atol=0.001
    )
    np.testing.assert_almost_equal(m_new.min_score, -13.6978, decimal=4)
    np.testing.assert_almost_equal(m_new.max_score, 4.1655, decimal=4)


def test_motif_repr(my_motif):
    assert "unnamed_motif_nACGTmskrwyn" == str(my_motif)


def test_information_content(pfm):
    m = Motif(pfm[:1])
    np.testing.assert_almost_equal(m.information_content, 0)
    m = Motif(pfm[:2])
    np.testing.assert_almost_equal(m.information_content, 2, decimal=2)
    m = Motif(pfm)
    np.testing.assert_almost_equal(m.information_content, 13.97, decimal=2)


def test_score_kmer(my_motif):
    with pytest.raises(ValueError):
        my_motif.score_kmer("ACG")

    np.testing.assert_almost_equal(
        my_motif.score_kmer("ATTTATTACCGT"), my_motif.min_score
    )
    np.testing.assert_almost_equal(
        my_motif.score_kmer("AACGTACGAACA"), my_motif.max_score
    )
    np.testing.assert_almost_equal(
        my_motif.score_kmer("CGTAGCGTTTAT"), -24.584, decimal=3
    )


def test_motif_shuffle(my_motif):
    shuffled = my_motif.shuffle()
    assert sorted(my_motif.consensus) == sorted(shuffled.consensus)

    check = [my_motif.shuffle().consensus != my_motif.consensus for _ in range(10)]
    assert np.any(check)


def test_motif_hash(my_motif):

    assert my_motif.hash == "fb543845f4e50436"

    # different id should not result in different hash
    my_motif.id = "Different ID"
    assert my_motif.hash == "fb543845f4e50436"

    # slight changes should not result in different hash
    ppm = my_motif.ppm
    for i in range(len(ppm)):
        ppm[i, np.random.randint(0, 3)] += np.random.random() / 1e4

    assert Motif(ppm=ppm).hash == "fb543845f4e50436"


@pytest.mark.parametrize("kind", ["information", "frequency", "ensembl", "energy"])
def test_plot_logo(my_motif, kind):
    for add_left in [0, 1]:
        my_motif.plot_logo(fname="test/test.png", kind=kind, add_left=add_left)
        assert os.path.exists("test/test.png")
        os.unlink("test/test.png")


def test_read_motifs_pfm(pfmfile2):
    with open(pfmfile2) as f:
        motifs = read_motifs(f, fmt="pwm")

    motif_ids = [m.id for m in motifs]
    assert 5 == len(motif_ids)
    assert [
        "M1500_1.01",
        "M5659_1.01",
        "M5669_1.01",
        "M5715_1.01",
        "M5717_1.01",
    ] == motif_ids


def test_read_motifs_jaspar(jaspar):
    with open(jaspar) as f:
        motifs = read_motifs(f, fmt="jaspar")

    my_motifs = ["MA0002.2", "MA0003.3", "MA0004.1", "MA0006.1"]

    my_lens = [6, 6, 11, 11]

    motif_ids = [m.id for m in motifs]
    assert 4 == len(motif_ids)
    assert my_motifs == motif_ids
    assert my_lens == sorted([len(m) for m in motifs])


def ppm_to_str(self):
    ppm = [[0.01, 0.01, 0.01, 0.97], [0.123, 0.456, 0.222, 0.199]]

    m = Motif(ppm)

    s2 = "0.01\t0.01\t0.01\t0.97\n0.12\t0.46\t0.22\t0.20"
    s3 = "0.010\t0.010\t0.010\t0.970\n0.123\t0.456\t0.222\t0.199"

    assert s2 == m._pwm_to_str(precision=2)
    assert s3 == m._pwm_to_str(precision=3)


def test_logodds_matrix():
    pwm = [[0.5, 0.4, 0.1, 0.0], [0.25, 0.25, 0.25, 0.25]]

    logodds = np.array(
        [
            [0.69813, 0.47623, -0.8916, -4.60517],
            [0.00995, 0.00995, 0.00995, 0.00995],
        ]
    )
    m = Motif(pwm)
    np.testing.assert_almost_equal(logodds, np.array(m.logodds), decimal=5)


def test_read_motifs(pfmfile2):

    # Read motifs from file
    motifs = read_motifs(pfmfile2, fmt="pwm")
    assert 5 == len(motifs)

    # Read motifs from file as dictionary
    motifs = read_motifs(pfmfile2, fmt="pwm", as_dict=True)
    assert 5 == len(motifs)
    assert type({}) == type(motifs)


def test_motif_export_import():
    pfm = [
        [120, 0, 0, 0],
        [120, 0, 0, 0],
        [0, 60, 60, 0],
        [0, 0, 0, 120],
        [0, 0, 0, 120],
    ]
    motif = Motif(pfm)
    motif.id = "test_motif"

    f = StringIO(motif.to_transfac())
    motif_from_file = read_motifs(f, fmt="transfac")[0]
    assert "AASTT" == motif_from_file.to_consensus().upper()
    assert "test_motif" == motif_from_file.id

    f = StringIO(motif.to_meme())
    motif_from_file = read_motifs(f, fmt="meme")[0]
    assert "AASTT" == motif_from_file.to_consensus().upper()
    assert "test_motif" == motif_from_file.id

    f = StringIO(motif.to_motevo())
    motif_from_file = read_motifs(f, fmt="transfac")[0]
    assert "AASTT" == motif_from_file.to_consensus().upper()
    assert "test_motif" == motif_from_file.id


def test_motif_from_alignment():
    align = "AACTT\n" "AAGTA\n" "AACTC\n" "AAGTG\n"
    f = StringIO(align)
    motif = read_motifs(f, fmt="align")[0]

    assert "AASTN" == motif.to_consensus().upper()


def test_read_motifs_xxmotifs():
    fname = "test/data/motifprogram/xxmotif.pwm"
    motifs = read_motifs(fname, fmt="xxmotif")

    assert 4 == len(motifs)
    assert 9 == len(motifs[-1])
    assert "RGGCAWGYC" == motifs[-1].to_consensus().upper()


def test6_pcc():
    pfm1 = [[5, 0, 0, 0], [0, 5, 0, 0], [0, 5, 0, 0], [0, 0, 0, 5]]
    pfm2 = [[5, 0, 0, 0], [0, 5, 0, 0], [0, 5, 0, 0], [0, 0, 0, 5]]

    m1 = Motif(pfm1)
    m2 = Motif(pfm2)

    assert 4 == m1.max_pcc(m2)[0]


def test_add_operator():
    m1 = motif_from_consensus("AAA")
    m2 = motif_from_consensus("TCG")

    assert (m1 + m2).consensus.upper() == "WMR"

    m2 = Motif([[0, 0, 0, 2], [0, 2, 0, 0], [0, 0, 2, 0]])
    assert (m1 + m2).consensus.upper() == "AAA"


def test_and_operator():
    m1 = motif_from_consensus("AAA")
    m2 = motif_from_consensus("TCG")

    assert (m1 & m2).consensus.upper() == "WMR"


def test_mul_operator():
    m1 = motif_from_consensus("AAA")
    m2 = motif_from_consensus("TCG")

    assert ((m1 * 10) + m2).consensus.upper() == "AAA"


def test_shift_operators():
    m = motif_from_consensus("AA")

    assert (m >> 1).consensus.upper() == "NAA"
    assert (m >> 2).consensus.upper() == "NNAA"
    assert (m << 1).consensus.upper() == "AAN"
    assert (m << 2).consensus.upper() == "AANN"


def test_invert_operator():
    m = motif_from_consensus("WMR")

    assert (~m).consensus.upper() == "YKW"
