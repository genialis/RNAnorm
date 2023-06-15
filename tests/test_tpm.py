import numpy as np
import pandas as pd
import pytest

from rnanorm import TPM
from rnanorm.annotation import GTF


@pytest.fixture
def expected(exp):
    return pd.DataFrame(
        [
            [83333.333, 83333.333, 83333.333, 166666.666, 583333.333],
            [83333.333, 83333.333, 83333.333, 166666.666, 583333.333],
            [45454.545, 45454.545, 45454.545, 90909.090, 772727.272],
            [142857.142, 142857.142, 142857.142, 285714.285, 285714.285],
        ],
        columns=exp.columns,
        index=exp.index,
    )


def test_tpm(exp, expected, gtf_file):
    transformer = TPM(gtf=gtf_file)

    # Numpy array is not allowed when using GTF
    with pytest.raises(ValueError, match=r"X should be a pandas.DataFrame"):
        exp_normalized = transformer.fit_transform(exp.values)

    # Default output is np.ndarray
    exp_normalized = transformer.fit_transform(exp)
    assert isinstance(exp_normalized, np.ndarray)
    np.testing.assert_array_almost_equal(
        exp_normalized,
        expected.to_numpy(),
        decimal=3,
    )

    # Set pandas output
    transformer.set_output(transform="pandas")
    exp_normalized = transformer.fit_transform(exp)
    assert isinstance(exp_normalized, pd.DataFrame)
    pd.testing.assert_frame_equal(exp_normalized, expected)

    # Validation refuses nan values
    exp.iloc[0, 0] = np.nan
    with pytest.raises(ValueError, match="Input X contains NaN."):
        transformer.fit_transform(exp)


def test_tpm_gene_lengths(exp, expected, gtf_file):
    gene_lengths = GTF(gtf=gtf_file).length

    # Numpy arrays are allowed when using gene lengths.
    # But numpy array must be subset ti have the same number of
    # genes as the gene lengths.
    transformer = TPM(gene_lengths=gene_lengths.to_numpy())
    exp_normalized = transformer.fit_transform(exp.values)
    np.testing.assert_array_almost_equal(
        exp_normalized,
        expected.to_numpy(),
        decimal=3,
    )

    # Work with pandas.
    transformer = TPM(gene_lengths=gene_lengths)
    transformer.set_output(transform="pandas")
    exp_normalized = transformer.fit_transform(exp)
    assert isinstance(exp_normalized, pd.DataFrame)
    pd.testing.assert_frame_equal(exp_normalized, expected)
