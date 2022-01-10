import numpy as np
import pandas as pd
import pytest

from rnanorm.normalization import fpkm


def test_fpkm_normalization():
    """Test FPKM formula implementation.

    FPKM expressions of a minimal example have been manually computed and
    compared against implementation.

    Expression data:

                        S1     S2
    ENSG00000136807   1000    400
    ENSG00000176903   2000  24000
    ENSG00000241490  17000    600

    Gene lengths:

    ENSG00000136807    1000
    ENSG00000176903    2000
    ENSG00000241490    5000

    Manually computed FPKM (rounded):

                         S1      S2
    ENSG00000136807       1       2
    ENSG00000176903       2     1.5
    ENSG00000241490       1       1
    """
    genes = ["ENSG00000136807", "ENSG00000176903", "ENSG00000241490"]
    gene_lengths = [1000, 2000, 5000]
    expressions = [[1000, 400], [2000, 24000], [17000, 600]]
    manually_computed_FPKM = [[50000, 16000], [50000, 480000], [170000, 4800]]

    X = pd.DataFrame(expressions, index=genes, columns=["S1", "S2"])
    y = pd.DataFrame(gene_lengths, index=genes, columns=["GENE_LENGTHS"])

    FPKM = fpkm(X, y)
    assert np.all(np.asarray(FPKM, dtype=int) == manually_computed_FPKM)  # Test Pandas array

    FPKM = fpkm(X.to_numpy(), y.to_numpy())
    assert np.all(FPKM.astype(int) == manually_computed_FPKM)  # Test Numpy array


def test_devision_by_zero():
    genes = ["ENSG00000136807", "ENSG00000176903", "ENSG00000241490"]
    gene_lengths = [1000, 2000, 3000]
    expressions = [[0], [0], [0]]
    expected_FPKM = [[0], [0], [0]]

    X = pd.DataFrame(expressions, index=genes, columns=["S1"])
    y = pd.DataFrame(gene_lengths, index=genes, columns=["GENE_LENGTHS"])

    FPKM = fpkm(X.to_numpy(), y.to_numpy())
    print(FPKM)

    assert np.all(FPKM.astype(int) == expected_FPKM)


def test_geneset_mismatch():
    genes = ["ENSG00000136807", "ENSG00000176903", "ENSG00000241490"]
    gene_lengths = [3000, 2000]
    expressions = [[0], [0], [0]]
    expected_FPKM = [[0], [0]]

    X = pd.DataFrame(expressions, index=genes, columns=["S1"])
    y = pd.DataFrame(gene_lengths, index=genes[:-1], columns=["GENE_LENGTHS"])

    with pytest.warns(RuntimeWarning, match="Geneset mismatch"):
        FPKM = fpkm(X, y)

    assert np.all(np.asarray(FPKM, dtype=int) == expected_FPKM)
