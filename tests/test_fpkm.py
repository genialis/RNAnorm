import numpy as np
import pandas as pd
import pytest

from rnanorm import FPKM


@pytest.fixture
def expected(exp):
    return pd.DataFrame(
        [
            [100000, 100000, 100000, 200000, np.nan],
            [100000, 100000, 100000, 200000, np.nan],
            [50000, 50000, 50000, 100000, np.nan],
            [200000, 200000, 200000, 400000, np.nan],
        ],
        columns=exp.columns,
        index=exp.index,
        dtype=np.float64,
    )


def test_fpkm(exp, expected, gtf_file):
    transformer = FPKM(gtf=gtf_file)

    # Default output is np.ndarray
    with pytest.warns(UserWarning, match=r"X contains .* genes that are not ."):
        exp_normalized = transformer.fit_transform(exp)
    assert isinstance(exp_normalized, np.ndarray)
    np.testing.assert_array_equal(exp_normalized, expected.to_numpy())

    # Check it possible to also get pandas output
    transformer.set_output(transform="pandas")
    with pytest.warns(UserWarning, match=r"X contains .* genes that are not ."):
        exp_normalized = transformer.fit_transform(exp)
    assert isinstance(exp_normalized, pd.DataFrame)
    pd.testing.assert_frame_equal(exp_normalized, expected)

    # Validation refuses nan values
    exp.iloc[0, 0] = np.nan
    with pytest.raises(ValueError, match="Input X contains NaN."):
        transformer.fit_transform(exp)
