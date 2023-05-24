import pathlib

import numpy as np
import pandas as pd

from rnanorm import TPM
from rnanorm.datasets import load_rnaseq_toy


def test_example():
    """Test that datasets.load_rnaseq_toy works as expected."""
    # Numpy array output
    ds = load_rnaseq_toy()
    assert isinstance(ds.exp, np.ndarray)
    assert isinstance(ds.gtf_path, pathlib.Path)

    # Numpy array output
    ds = load_rnaseq_toy(as_frame=True)
    assert isinstance(ds.exp, pd.DataFrame)

    result = TPM(ds.gtf_path).fit_transform(ds.exp)

    assert isinstance(result, np.ndarray)
    assert result.shape == (4, 5)
