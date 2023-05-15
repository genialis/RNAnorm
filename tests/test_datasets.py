import numpy as np

from rnanorm import TPM
from rnanorm.datasets import load_example


def test_example():
    """Test that datasets.load_example works as expected."""
    example = load_example()
    result = TPM(example.gtf_path).fit_transform(example.exp)

    assert isinstance(result, np.ndarray)
    assert result.shape == (4, 5)
