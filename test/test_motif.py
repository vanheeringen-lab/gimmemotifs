import pytest
import numpy as np

from gimmemotifs.motif import Motif

@pytest.fixture
def pfm():
    return  [
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
    np.testing.assert_allclose(m_new.ppm, [[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], atol=0.001)
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

    np.testing.assert_almost_equal(my_motif.score_kmer("ATTTATTACCGT"), my_motif.min_score)
    np.testing.assert_almost_equal(my_motif.score_kmer("AACGTACGAACA"), my_motif.max_score)
    np.testing.assert_almost_equal(my_motif.score_kmer("CGTAGCGTTTAT"), -24.584, decimal=3)


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

    assert Motif(ppm).hash == "fb543845f4e50436"



# class TestMotif(unittest.TestCase):
#     """ A test class for Motif """

#     def setUp(self):
#         self.data_dir = "test/data/motif"
#         self.pwm = os.path.join(self.data_dir, "test.pwm")
#         self.pwm2 = "test/data/pwms/motifs.pwm"
#         self.jaspar = "test/data/pwms/test.jaspar"
#         self.pfm = [
#             [1, 0, 0, 0],
#             [0, 1, 0, 0],
#             [0, 0, 1, 0],
#             [0, 0, 0, 1],
#             [1, 1, 0, 0],
#             [0, 1, 1, 0],
#             [0, 0, 1, 1],
#             [2, 0, 2, 0],
#             [3, 0, 0, 3],
#             [0, 1, 0, 1],
#         ]


#     def test5_motif_to_img(self):
#         """ Motif to img """
#         seqlogo = which("seqlogo")
#         if seqlogo:
#             m = Motif(self.pfm)
#             m.to_img("test/test.png", fmt="png", seqlogo=seqlogo)
#             self.assertTrue(os.path.exists("test/test.png"))
#             os.unlink("test/test.png")
#         else:
#             print("seqlogo not found, skipping.")

#     def test6_pcc(self):
#         pfm1 = [[5, 0, 0, 0], [0, 5, 0, 0], [0, 5, 0, 0], [0, 0, 0, 5]]
#         pfm2 = [[5, 0, 0, 0], [0, 5, 0, 0], [0, 5, 0, 0], [0, 0, 0, 5]]

#         m1 = Motif(pfm1)
#         m2 = Motif(pfm2)

#         self.assertEqual(4, m1.max_pcc(m2)[0])

#     def test7__read_motifs_pwm(self):
#         with open(self.pwm2) as f:
#             motifs = read_motifs(f, fmt="pwm")

#         motif_ids = [m.id for m in motifs]
#         self.assertEqual(5, len(motif_ids))
#         self.assertEqual(
#             ["M1500_1.01", "M5659_1.01", "M5669_1.01", "M5715_1.01", "M5717_1.01"],
#             motif_ids,
#         )

#     def test7__read_motifs_jaspar(self):
#         with open(self.jaspar) as f:
#             motifs = read_motifs(f, fmt="jaspar")

#         my_motifs = ["MA0002.2", "MA0003.3", "MA0004.1", "MA0006.1"]

#         my_lens = [6, 6, 11, 11]

#         motif_ids = [m.id for m in motifs]
#         self.assertEqual(4, len(motif_ids))
#         self.assertEqual(my_motifs, motif_ids)
#         self.assertEqual(my_lens, sorted([len(m) for m in motifs]))

#     def test8_pwm_to_str(self):
#         pwm = [[0.01, 0.01, 0.01, 0.97], [0.123, 0.456, 0.222, 0.199]]

#         m = Motif(pwm)

#         s2 = "0.01\t0.01\t0.01\t0.97\n0.12\t0.46\t0.22\t0.20"
#         s3 = "0.010\t0.010\t0.010\t0.970\n0.123\t0.456\t0.222\t0.199"

#         self.assertEqual(s2, m._pwm_to_str(precision=2))
#         self.assertEqual(s3, m._pwm_to_str(precision=3))

#     def test8_pwm_to_str_hash(self):
#         pwm = [[0.01, 0.01, 0.01, 0.97], [0.123, 0.456, 0.222, 0.199]]
#         m = Motif(pwm)
#         h = "1f260320cac8c26a"
#         self.assertEqual(h, m.hash())

#         pwm = [
#             [0.010000, 0.010000, 0.010000, 0.970000],
#             [0.12300, 0.45600, 0.22200, 0.19900],
#         ]
#         m = Motif(pwm)
#         self.assertEqual(h, m.hash())

#     def test9_logodds_matrix(self):
#         pwm = [[0.5, 0.4, 0.1, 0.0], [0.25, 0.25, 0.25, 0.25]]

#         logodds = np.array(
#             [
#                 [0.69813, 0.47623, -0.89160, -4.60517],
#                 [0.00995, 0.00995, 0.00995, 0.00995],
#             ]
#         )
#         m = Motif(pwm)
#         np.testing.assert_almost_equal(logodds, np.array(m.logodds), decimal=5)

#     def test10_read_motifs(self):

#         # Read motifs from file
#         motifs = read_motifs(self.pwm2, fmt="pwm")
#         self.assertEqual(5, len(motifs))

#         # Read motifs from file as dictionary
#         motifs = read_motifs(self.pwm2, fmt="pwm", as_dict=True)
#         self.assertEqual(5, len(motifs))
#         self.assertEqual(type({}), type(motifs))


#     def test_motif_export_import(self):
#         pfm = [
#             [120, 0, 0, 0],
#             [120, 0, 0, 0],
#             [0, 60, 60, 0],
#             [0, 0, 0, 120],
#             [0, 0, 0, 120],
#         ]
#         motif = Motif(pfm)
#         motif.id = "test_motif"

#         f = StringIO(motif.to_transfac())
#         motif_from_file = read_motifs(f, fmt="transfac")[0]
#         self.assertEqual("AASTT", motif_from_file.to_consensus().upper())
#         self.assertEqual("test_motif", motif_from_file.id)

#         f = StringIO(motif.to_meme())
#         motif_from_file = read_motifs(f, fmt="meme")[0]
#         self.assertEqual("AASTT", motif_from_file.to_consensus().upper())
#         self.assertEqual("test_motif", motif_from_file.id)

#         f = StringIO(motif.to_motevo())
#         motif_from_file = read_motifs(f, fmt="transfac")[0]
#         self.assertEqual("AASTT", motif_from_file.to_consensus().upper())
#         self.assertEqual("test_motif", motif_from_file.id)

#     def test_motif_from_alignment(self):
#         align = "AACTT\n" "AAGTA\n" "AACTC\n" "AAGTG\n"
#         f = StringIO(align)
#         motif = read_motifs(f, fmt="align")[0]

#         self.assertEqual("AASTN", motif.to_consensus().upper())

#     def test_read_motifs_xxmotifs(self):
#         fname = "test/data/motifprogram/xxmotif.pwm"
#         motifs = read_motifs(fname, fmt="xxmotif")

#         self.assertEqual(4, len(motifs))
#         self.assertEqual(9, len(motifs[-1]))
#         self.assertEqual("RGGCAWGYC", motifs[-1].to_consensus().upper())

#     def tearDown(self):
#         pass


# if __name__ == "__main__":
#     unittest.main()
