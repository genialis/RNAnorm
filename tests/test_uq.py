import numpy as np

from rnanorm import UQ


def test_uq(between_sample_data):
    """Test UQ normalization."""
    data, expected_factors, expected_data, _ = between_sample_data

    # Simple case: Fit some data, transform same data.
    transformer = UQ()
    transformer.fit(data)
    factors = transformer.get_norm_factors(data)
    transformed_data = transformer.transform(data)

    np.testing.assert_allclose(factors, expected_factors, rtol=1e-4)
    np.testing.assert_allclose(transformed_data, expected_data, rtol=1e-4)

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
