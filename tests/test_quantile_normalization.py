import numpy as np
import pandas as pd

from rnanorm.normalization import quantile


def test_quantile_normalization():
    """Test quantile normalization implementation.

    Quantile normalization of a minimal example has been computed manually and
    compared against implementation.

    Expression data:

        S1   S2   S3
    A    5    4    3
    B    2    1    4
    C    3    4    6
    D    4    2    8

    Manually computed quantile normalization (rounded):

           S1      S2      S3
    A    5.67    5.17    2.00
    B    2.00    2.00    3.00
    C    3.00    5.17    4.67
    D    4.67    3.00    5.67
    """
    expressions = [[5, 4, 3], [2, 1, 4], [3, 4, 6], [4, 2, 8]]
    manually_computed_qt = [[5.67, 5.17, 2.00], [2.00, 2.00, 3.00], [3.00, 5.17, 4.67], [4.67, 3.00, 5.67]]

    X = pd.DataFrame(expressions)

    # Test Pandas DataFrame
    qtX = quantile(X)
    qtX_rounded = np.round(np.asarray(qtX, dtype=float), 2)
    assert np.all(qtX_rounded == manually_computed_qt)
    assert type(qtX) == pd.DataFrame  # Check DataFrame returned for DataFrame input

    # Test Numpy array
    qtX = quantile(X.to_numpy())
    qtX_rounded = np.round(qtX, 2)
    assert np.all(qtX_rounded == manually_computed_qt)
    assert type(qtX) == np.ndarray  # Check Numpy array returned for Numpy array input

    # Test Python array
    qtX = quantile(list(X.to_numpy()))
    qtX_rounded = np.round(qtX, 2)
    assert np.all(qtX_rounded == manually_computed_qt)  # Test Numpy array
    assert type(qtX) == np.ndarray  # Check Numpy array returned for Python array input
