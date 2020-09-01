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


def test_quantile_normalization_special():
    """
    Test quantile normalization for special cases.

    Check if some special cases are handled the same as in the R implementation.
    """
    # tie of more then two numbers
    x = [[5, 4, 3], [2, 4, 4], [3, 4, 7], [4, 2, 8]]
    x_man = [[5.67, 5.00, 2.33], [2.33, 5.00, 3.67], [3.67, 5.00, 5.00], [5.00, 2.33, 5.67]]
    x_qt = np.round(quantile(x), 2)
    assert np.all(x_qt == x_man)

    # all zeros
    x = [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]]
    x_man = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    x_qt = np.round(quantile(x), 2)
    assert np.all(x_qt == x_man)

    # rank is all zeros
    x = [[5, 4, 0], [0, 0, 4], [3, 4, 6], [4, 2, 8]]
    x_man = [[5.67, 5.17, 0.00], [0.00, 0.00, 3.00], [3.00, 5.17, 4.67], [4.67, 3.00, 5.67]]
    x_qt = np.round(quantile(x), 2)
    assert np.all(x_qt == x_man)

    # column is all zeros
    x = [[0, 4, 3], [0, 1, 4], [0, 4, 6], [0, 2, 8]]
    x_man = [[2.67, 3.67, 1.33], [2.67, 1.33, 2.00], [2.67, 3.67, 3.33], [2.67, 2.00, 4.00]]
    x_qt = np.round(quantile(x), 2)
    assert np.all(x_qt == x_man)
