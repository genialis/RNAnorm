import numpy as np
import pandas as pd

from rnanorm.normalization import tpm


def test_tpm_normalization():
    """Test TPM formula implementation.

    TPM expressions of a minimal example have been manually computed and
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

    Manually computed TPM (rounded):

                         S1      S2
    ENSG00000136807   20704   17400
    ENSG00000176903  959266  957349
    ENSG00000241490   20029   25250
    """
    genes = ["ENSG00000136807", "ENSG00000176903", "ENSG00000241490"]
    gene_lengths = [3000, 2590, 3101]
    expressions = [[500, 400], [20000, 19000], [500, 600]]
    manually_computed_TPM = [[20704, 17400], [959266, 957349], [20029, 25250]]

    X = pd.DataFrame(expressions, index=genes, columns=["S1", "S2"])
    y = pd.DataFrame(gene_lengths, index=genes, columns=["GENE_LENGTHS"])

    TPM = tpm(X, y)
    assert np.all(np.asarray(TPM, dtype=int) == manually_computed_TPM)  # Test Pandas array

    TPM = tpm(X.to_numpy(), y.to_numpy())
    assert np.all(TPM.astype(int) == manually_computed_TPM)  # Test Numpy array


def test_devision_by_zero():
    genes = ["ENSG00000136807", "ENSG00000176903", "ENSG00000241490"]
    gene_lengths = [3000, 2590, 3101]
    expressions = [[0], [0], [0]]
    expected_TPM = [[0], [0], [0]]

    X = pd.DataFrame(expressions, index=genes, columns=["S1"])
    y = pd.DataFrame(gene_lengths, index=genes, columns=["GENE_LENGTHS"])

    TPM = tpm(X.to_numpy(), y.to_numpy())

    assert np.all(TPM.astype(int) == expected_TPM)
