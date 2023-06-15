import numpy as np
import pandas as pd
import pytest

from rnanorm import CTF


@pytest.fixture
def expected_ctf(exp):
    return pd.DataFrame(
        [
            [200, 300, 500, 2000, 7000],
            [400, 600, 1000, 4000, 14000],
            [400, 600, 1000, 4000, 34000],
            [100, 150, 250, 1000, 1000],
        ],
        index=exp.index,
        columns=exp.columns,
        dtype=np.float64,
    )


def test_ctf(exp, expected_factors, expected_ctf):
    """Test UQ normalization."""
    # Simple case: Fit some data, transform same data.
    transformer = CTF()
    transformer.fit(exp)
    factors = transformer.get_norm_factors(exp)
    transformed_data = transformer.transform(exp)

    np.testing.assert_allclose(factors, expected_factors, rtol=1e-4)
    np.testing.assert_allclose(transformed_data, expected_ctf, rtol=1e-4)

    # Advanced case: Fit some samples, transform different ones.
    samples_to_fit = ["Sample_1", "Sample_3", "Sample_4"]
    exp1 = exp.loc[samples_to_fit]
    transformer.fit(exp1)
    factors = transformer.get_norm_factors(exp1)
    transformed_data = transformer.transform(exp.loc[["Sample_2"]])

    np.testing.assert_allclose(factors, [1, 0.5, 2.0], rtol=1e-4)
    np.testing.assert_allclose(
        transformed_data,
        expected_ctf.loc[["Sample_2"]],
        rtol=1e-3,
    )
