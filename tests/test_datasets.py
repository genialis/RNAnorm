import pathlib

import numpy as np
import pandas as pd

from rnanorm import TPM
from rnanorm.datasets import load_gtex, load_rnaseq_toy


def test_rnaseq_toy():
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


def test_gtex():
    """Test that datasets.load_gtex works as expected."""
    # Numpy array output
    ds = load_gtex()
    assert isinstance(ds.exp, np.ndarray)
    assert isinstance(ds.gtf_path, pathlib.Path)

    # Pandas array output
    ds = load_gtex(as_frame=True)
    assert isinstance(ds.exp, pd.DataFrame)
    assert ds.exp.shape == (30, 818)
    assert ds.exp.index[0] == "GTEX-111CU-0326-SM-5GZXO"
    assert ds.exp.columns[0] == "ENSG00000141956.13"

    result = TPM(ds.gtf_path).set_output(transform="pandas").fit_transform(ds.exp)

    assert isinstance(result, pd.DataFrame)
    assert result.shape == (30, 818)
