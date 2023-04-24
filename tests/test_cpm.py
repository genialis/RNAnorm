import numpy as np
import pandas as pd
import pytest

from rnanorm import CPM


def test_cpm():
    exp = pd.DataFrame(
        [
            [5, 5, 0],
            [4, 1, 0],
        ],
        index=["SAMPLE_1", "SAMPLE_2"],
        columns=["GENE_2", "GENE_1", "GENE_3"],
    )
    expected = pd.DataFrame(
        [
            [500000.0, 500000.0, 0],
            [800000.0, 200000.0, 0.0000],
        ],
        index=exp.index,
        columns=exp.columns,
    )
    transformer = CPM()

    # Default output is np.ndarray
    exp_normalized = transformer.fit_transform(exp)
    assert isinstance(exp_normalized, np.ndarray)
    np.testing.assert_array_equal(exp_normalized, expected.to_numpy())

    # Test nested Python list
    exp_normalized = transformer.fit_transform(exp.values.tolist())
    assert isinstance(exp_normalized, np.ndarray)
    assert exp_normalized.tolist() == expected.values.tolist()

    # Check it possible to also get pandas output
    transformer.set_output(transform="pandas")
    exp_normalized = transformer.fit_transform(exp)
    assert isinstance(exp_normalized, pd.DataFrame)
    pd.testing.assert_frame_equal(exp_normalized, expected)

    # Validation refuses nan values
    exp.iloc[0, 0] = np.nan
    with pytest.raises(ValueError, match="Input X contains NaN."):
        transformer.fit_transform(exp)
