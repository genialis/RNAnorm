import numpy as np

from rnanorm import CUF


def test_cuf(between_sample_data):
    """Test CUF normalization."""
    data, expected_factors, _, expected_data = between_sample_data

    transformer = CUF()
    transformer.fit(data)
    factors = transformer.get_norm_factors(data)
    transformed_data = transformer.transform(data)

    np.testing.assert_allclose(factors, expected_factors, rtol=1e-3)
    np.testing.assert_allclose(transformed_data, expected_data, rtol=1e-4)

    # Advanced case: Fit some samples, transform different ones.
    transformer.fit(data.iloc[1:])
    factors = transformer.get_norm_factors(data.iloc[:1])
    transformed_data = transformer.transform(data.iloc[:1])

    np.testing.assert_allclose(factors, [1.4142], rtol=1e-4)
    np.testing.assert_allclose(
        transformed_data,
        np.array([[0.0, 0.7071, 7.071, 70.71, 707.1, 7071]]),
        rtol=1e-3,
    )
