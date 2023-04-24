import numpy as np
import pandas as pd
import pytest

from rnanorm import LibrarySize


def test_lib_size():
    exp = pd.DataFrame(
        [
            [5, 5, 0],
            [0, 0, 0],
            [4, 1, 0],
        ],
        index=["SAMPLE_1", "SAMPLE_2", "SAMPLE_3"],
        columns=["GENE_1", "GENE_2", "GENE_3"],
        dtype=float,
    )
    expected = pd.DataFrame(
        [10.0, 0.0, 5.0],
        index=exp.index,
        columns=["Library size"],
    )

    transformer = LibrarySize()

    # Default output is np.ndarray
    lib_size = transformer.fit_transform(exp)
    assert isinstance(lib_size, np.ndarray)
    np.testing.assert_array_equal(lib_size, expected.to_numpy()[:, 0])

    # Check it possible to also get pandas output
    transformer.set_output(transform="pandas")
    lib_size = transformer.fit_transform(exp)
    assert isinstance(lib_size, pd.DataFrame)
    pd.testing.assert_frame_equal(lib_size, expected)

    # Validation refuses nan values
    exp.iloc[0, 0] = np.nan
    with pytest.raises(ValueError, match="Input X contains NaN."):
        transformer.fit_transform(exp)
