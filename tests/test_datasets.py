import pandas as pd

from rnanorm import TPM
from rnanorm.datasets import load_gtex, load_toy_data


def test_toy_data():
    """Test that datasets.load_toy_data works as expected."""
    ds = load_toy_data()
    assert isinstance(ds.exp, pd.DataFrame)
    assert ds.exp.shape == (4, 5)
    assert ds.exp.index[0] == "Sample_1"
    assert ds.exp.columns[0] == "Gene_1"

    result = TPM(ds.gtf_path).set_output(transform="pandas").fit_transform(ds.exp)

    assert isinstance(result, pd.DataFrame)
    assert result.shape == (4, 5)


def test_gtex():
    """Test that datasets.load_gtex works as expected."""
    ds = load_gtex()
    assert isinstance(ds.exp, pd.DataFrame)
    assert ds.exp.shape == (30, 818)
    assert ds.exp.index[0] == "GTEX-111CU-0326-SM-5GZXO"
    assert ds.exp.columns[0] == "ENSG00000141956.13"

    result = TPM(ds.gtf_path).set_output(transform="pandas").fit_transform(ds.exp)

    assert isinstance(result, pd.DataFrame)
    assert result.shape == (30, 818)
