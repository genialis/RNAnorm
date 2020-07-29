import numpy as np
import pandas as pd

from rnanorm.normalization import cpm


def test_cpm_normalization():
    """Test CPM formula implementation.

    CPM expressions of a minimal example have been manually computed and
    compared against implementation.

    Expression data:

                        S1     S2
    ENSG00000136807    500    400
    ENSG00000176903  20000  19000
    ENSG00000241490    500    600

    Gene lengths:

    ENSG00000136807  3000
    ENSG00000176903  2590
    ENSG00000241490  3101

    Manually computed CPM (rounded):

                         S1      S2
    ENSG00000136807   23809   20000
    ENSG00000176903  952380  950000
    ENSG00000241490   23809   30000
    """
    expressions = [[500, 400], [20000, 19000], [500, 600]]
    manually_computed_CPM = [[23809, 20000], [952380, 950000], [23809, 30000]]

    X = pd.DataFrame(expressions, columns=["S1", "S2"])

    # Test Pandas DataFrame
    CPM = cpm(X)
    assert np.all(np.asarray(CPM, dtype=int) == manually_computed_CPM)
    assert type(CPM) == pd.DataFrame  # Check DataFrame returned for DataFrame input

    # Test Numpy array
    CPM = cpm(X.to_numpy())
    assert np.all(CPM.astype(int) == manually_computed_CPM)
    assert type(CPM) == np.ndarray  # Check Numpy array returned for Numpy array input

    # Test Python array
    CPM = cpm(list(X.to_numpy()))
    assert np.all(CPM.astype(int) == manually_computed_CPM)  # Test Numpy array
    assert type(CPM) == np.ndarray  # Check Numpy array returned for Python array input


def test_devision_by_zero():
    expressions = [[0], [0], [0]]
    expected_CPM = [[0], [0], [0]]

    X = pd.DataFrame(expressions, columns=["S1"])

    CPM = cpm(X)

    assert np.all(np.asarray(CPM, dtype=int) == expected_CPM)
