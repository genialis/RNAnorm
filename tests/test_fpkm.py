import numpy as np
import pandas as pd
import pytest

from rnanorm import FPKM
from rnanorm.annotation import GTF


@pytest.fixture
def expected(exp):
    return pd.DataFrame(
        [
            [100000, 100000, 100000, 200000, 700000],
            [100000, 100000, 100000, 200000, 700000],
            [50000, 50000, 50000, 100000, 850000],
            [200000, 200000, 200000, 400000, 400000],
        ],
        columns=exp.columns,
        index=exp.index,
        dtype=np.float64,
    )


def test_fpkm(exp, expected, gtf_file):
    transformer = FPKM(gtf=gtf_file)

    # Numpy array is not allowed when using GTF
    with pytest.raises(ValueError, match=r"X should be a pandas.DataFrame"):
        exp_normalized = transformer.fit_transform(exp.values)

    # Default output is np.ndarray
    exp_normalized = transformer.fit_transform(exp)
    assert isinstance(exp_normalized, np.ndarray)
    np.testing.assert_array_equal(exp_normalized, expected.to_numpy())

    # Set pandas output
    transformer.set_output(transform="pandas")
    exp_normalized = transformer.fit_transform(exp)
    assert isinstance(exp_normalized, pd.DataFrame)
    pd.testing.assert_frame_equal(exp_normalized, expected)

    # Validation refuses nan values
    exp.iloc[0, 0] = np.nan
    with pytest.raises(ValueError, match="Input X contains NaN."):
        transformer.fit_transform(exp)


def test_fpkm_gene_lengths(exp, expected, gtf_file):
    gene_lengths = GTF(gtf=gtf_file).length

    # Numpy arrays are allowed when using gene lengths.
    # But numpy array must be subset ti have the same number of
    # genes as the gene lengths.
    transformer = FPKM(gene_lengths=gene_lengths.to_numpy())
    exp_normalized = transformer.fit_transform(exp.values)
    np.testing.assert_array_equal(exp_normalized, expected.to_numpy())

    # Work with pandas.
    transformer = FPKM(gene_lengths=gene_lengths)
    transformer.set_output(transform="pandas")
    exp_normalized = transformer.fit_transform(exp)
    assert isinstance(exp_normalized, pd.DataFrame)
    pd.testing.assert_frame_equal(exp_normalized, expected)
