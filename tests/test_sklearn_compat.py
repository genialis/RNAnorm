import numpy as np
import pandas as pd
from sklearn import config_context
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from rnanorm import CPM, CTF, CUF, FPKM, TMM, TPM, UQ
from rnanorm.datasets import load_toy_data


def test_grid_search():
    """Test compatibility of all methods with sklearn machinery."""
    ds = load_toy_data()
    X = ds.exp
    y = pd.Series([0, 0, 1, 1], index=X.index)
    pipeline = Pipeline(
        steps=[
            ("normalization", CPM()),
            ("scaler", StandardScaler()),
            ("classifier", LogisticRegression()),
        ]
    )
    params = {
        "normalization": [
            CPM(),
            FPKM(gtf=ds.gtf_path),
            TPM(gtf=ds.gtf_path),
            UQ(),
            CUF(),
            TMM(),
            CTF(),
        ],
    }
    search = GridSearchCV(pipeline, params, cv=2, refit=False)
    search.fit(X, y)
    results = pd.DataFrame(search.cv_results_)
    assert results.shape[0] == 7


def test_set_output():
    """Ensure set_output behaves as expected."""
    ds = load_toy_data()

    method__has_factors = (
        (CPM(), False),
        (FPKM(gtf=ds.gtf_path), False),
        (TPM(gtf=ds.gtf_path), False),
        (UQ(), True),
        (CUF(), True),
        (TMM(), True),
        (CTF(), True),
    )

    for method, has_factors in method__has_factors:
        # No configuration returns np.arrays
        tmm = TMM()
        result = tmm.fit_transform(ds.exp)
        factors = tmm.get_norm_factors(ds.exp)
        assert isinstance(result, np.ndarray)
        assert isinstance(factors, np.ndarray)

        # Explicit global config should return corresponding objects
        with config_context(transform_output="default"):
            result = tmm.fit_transform(ds.exp)
            factors = tmm.get_norm_factors(ds.exp)
            assert isinstance(result, np.ndarray)
            assert isinstance(factors, np.ndarray)

        with config_context(transform_output="pandas"):
            result = tmm.fit_transform(ds.exp)
            factors = tmm.get_norm_factors(ds.exp)
            assert isinstance(result, pd.DataFrame)
            assert isinstance(factors, pd.Series)

        # Explicit pandas config should return pandas objects
        method.set_output(transform="pandas")
        result = method.fit_transform(ds.exp)
        assert isinstance(result, pd.DataFrame)
        if has_factors:
            factors = method.get_norm_factors(ds.exp)
            assert isinstance(factors, pd.Series)
