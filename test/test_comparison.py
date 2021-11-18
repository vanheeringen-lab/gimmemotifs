import numpy as np
import pytest

from gimmemotifs import comparison


@pytest.fixture
def pfm1():
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
def pfm2():
    return [
        [12, 0, 0, 0],
        [0, 12, 0, 0],
        [0, 0, 12, 0],
        [0, 0, 0, 12],
    ]


def test_make_equal_length(pfm1, pfm2):
    r1, r2 = comparison.make_equal_length(pfm1, pfm2, pos=1)

    assert len(r1) == len(r2)
    np.testing.assert_allclose(r2[0], [0.25, 0.25, 0.25, 0.25])
    np.testing.assert_allclose(r2[-1], [0.25, 0.25, 0.25, 0.25])

    r1, r2 = comparison.make_equal_length(np.array(pfm1), np.array(pfm2), pos=1)

    assert len(r1) == len(r2)
    np.testing.assert_allclose(r2[0], [0.25, 0.25, 0.25, 0.25])
    np.testing.assert_allclose(r2[-1], [0.25, 0.25, 0.25, 0.25])


def test_motif_comparer():
    mc = comparison.MotifComparer()

    pwm = "test/data/pwmscan/TATA.pwm"
    ret = mc.get_closest_match(pwm, parallel=True)

    assert "TATA-box" in ret

    match = ret["TATA-box"]
    assert "GM.5.0.TBP.0001" == match[0]

    scores = match[1]
    np.testing.assert_almost_equal(-0.1008, scores[0], 4)
    assert 0 == scores[1]
    assert 1 == scores[2]
    np.testing.assert_almost_equal(3.1666e-8, scores[3])
