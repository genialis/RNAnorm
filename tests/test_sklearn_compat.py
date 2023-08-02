import pandas as pd
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
