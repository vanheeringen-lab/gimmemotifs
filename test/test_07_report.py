import numpy as np
import pandas as pd

from gimmemotifs.report import ExtraStyler


def test_extrastyler():
    # this test should catch issues with pandas Styler methods
    # non_reducing_slice
    # _translate()
    df = pd.DataFrame(np.random.randn(4, 2), columns=["a", "b"])
    ExtraStyler(df).to_html(
        font="Nunito Sans",
        col_heading_style={
            "name": "col_heading",
            "props": [("border-bottom", "1px solid #e0e0e0")],
        },
    )
