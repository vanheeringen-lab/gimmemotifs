import pandas as pd
import numpy as np
from gimmemotifs.report import ExtraStyler


def test_extrastyler():
    # this test should catch issues with pandas Styler methods
    # non_reducing_slice
    # _translate()
    df = pd.DataFrame(np.random.randn(4, 2), columns=["a", "b"])
    ExtraStyler(df).render()
