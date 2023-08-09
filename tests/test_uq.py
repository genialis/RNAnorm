from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from rnanorm import UQ
from rnanorm.datasets import load_gtex


@pytest.fixture
def expected_uq(exp):
    return pd.DataFrame(
        [
            [20000, 30000, 50000, 200000, 700000],
            [20000, 30000, 50000, 200000, 700000],
            [20000, 30000, 50000, 200000, 1700000],
            [20000, 30000, 50000, 200000, 200000],
        ],
        index=exp.index,
        columns=exp.columns,
        dtype=np.float64,
    )


def test_uq(exp, expected_factors, expected_uq):
    """Test UQ normalization."""
    # Simple case: Fit some data, transform same data.
    transformer = UQ()
    transformer.fit(exp)
    factors = transformer.get_norm_factors(exp)
    transformed_data = transformer.transform(exp)

    np.testing.assert_allclose(factors, expected_factors, rtol=1e-4)
    np.testing.assert_allclose(transformed_data, expected_uq, rtol=1e-4)

    # Advanced case: Fit some samples, transform different ones.
    samples_to_fit = ["Sample_1", "Sample_3", "Sample_4"]
    exp1 = exp.loc[samples_to_fit]
    transformer.fit(exp1)
    factors = transformer.get_norm_factors(exp1)
    transformed_data = transformer.transform(exp.loc[["Sample_2"]])

    np.testing.assert_allclose(factors, [1, 0.5, 2.0], rtol=1e-4)
    np.testing.assert_allclose(
        transformed_data,
        expected_uq.loc[["Sample_2"]],
        rtol=1e-3,
    )


def test_uq_rnanorm_edger():
    """Test our results against EdgeR's UQ implementation."""
    # Load saved scaling factors from edgeR
    files_dir = Path(__file__).parent / "files"
    edger_factors = pd.read_csv(
        files_dir / "gtex.lung.30.chr21.uq_factors_edgeR.tsv",
        sep="\t",
        comment="#",
        index_col=0,
        usecols=[0, 2, 3],
    )

    # Compute scaling factors here
    ds = load_gtex()
    rnanorm_factors = UQ().fit(ds.exp).get_norm_factors(ds.exp)

    # Compare loaded and computed scaling factors
    np.testing.assert_array_almost_equal(
        edger_factors["norm.factors"].values,
        rnanorm_factors,
        decimal=14,
    )
