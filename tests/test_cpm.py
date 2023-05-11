import numpy as np
import pandas as pd
import pytest

from rnanorm import CPM


@pytest.fixture
def expected_cpm(exp):
    return pd.DataFrame(
        [
            [20000, 30000, 50000, 200000, 700000],
            [20000, 30000, 50000, 200000, 700000],
            [10000, 15000, 25000, 100000, 850000],
            [40000, 60000, 100000, 400000, 400000],
        ],
        index=exp.index,
        columns=exp.columns,
        dtype=np.float64,
    )


def test_cpm(exp, expected_cpm):
    transformer = CPM()

    # Default output is np.ndarray
    exp_normalized = transformer.fit_transform(exp)
    assert isinstance(exp_normalized, np.ndarray)
    np.testing.assert_array_equal(exp_normalized, expected_cpm.to_numpy())

    # Test nested Python list
    exp_normalized = transformer.fit_transform(exp.values.tolist())
    assert isinstance(exp_normalized, np.ndarray)
    assert exp_normalized.tolist() == expected_cpm.values.tolist()

    # Check it possible to also get pandas output
    transformer.set_output(transform="pandas")
    exp_normalized = transformer.fit_transform(exp)
    assert isinstance(exp_normalized, pd.DataFrame)
    pd.testing.assert_frame_equal(exp_normalized, expected_cpm)

    # Validation refuses nan values
    exp.iloc[0, 0] = np.nan
    with pytest.raises(ValueError, match="Input X contains NaN."):
        transformer.fit_transform(exp)
