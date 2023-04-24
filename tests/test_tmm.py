import numpy as np
import pytest

from rnanorm import TMM


def test_tmm(between_sample_data):
    """Test TMM normalization."""
    data, expected_factors, expected_data, _ = between_sample_data

    transformer = TMM()
    transformer.fit(data)
    factors = transformer.get_norm_factors(data)
    transformed_data = transformer.transform(data)

    assert data.shape == transformed_data.shape
    np.testing.assert_allclose(transformed_data, expected_data, rtol=1e-4)
    # Geometric mean of factors is 1:
    assert np.exp(np.mean(np.log(factors))) == pytest.approx(1)
    np.testing.assert_allclose(factors, expected_factors, rtol=1e-4)

    # Advanced case: Fit some samples, transform different ones.
    transformer.fit(data.iloc[1:])
    factors = transformer.get_norm_factors(data.iloc[:1])
    transformed_data = transformer.transform(data.iloc[:1])

    np.testing.assert_allclose(factors, [1.4142], rtol=1e-4)
    np.testing.assert_allclose(
        transformed_data,
        np.array([[0.0, 63.64, 636.4, 6364, 63640, 636400]]),
        rtol=1e-3,
    )
