import numpy as np
import pandas as pd
import pytest

from rnanorm import TPM


@pytest.mark.skip(reason="Skip - gtf_file fixture is synced with test data.")
def test_tpm(gtf_file):
    exp = pd.DataFrame(
        [
            [200000, 600000, 1.0],
            [400000, 200000, 0],
            [1500000, 1000000, 0],
        ],
        index=["SAMPLE_1", "SAMPLE_2", "SAMPLE_3"],
        columns=["GENE_1", "GENE_2", "GENE_3"],
    )

    expected = pd.DataFrame(
        [
            [400000.0, 600000.0, np.nan],
            [800000.0, 200000.0, np.nan],
            [750000.0, 250000.0, np.nan],
        ],
        index=exp.index,
        columns=exp.columns,
    )

    transformer = TPM(gtf=gtf_file)

    # Default output is np.ndarray
    with pytest.warns(UserWarning, match=r"X contains .* genes that are not ."):
        exp_normalized = transformer.fit_transform(exp)
    assert isinstance(exp_normalized, np.ndarray)
    np.testing.assert_array_equal(exp_normalized, expected.to_numpy())

    # Check it possible to also get pandas output
    transformer.set_output(transform="pandas")
    with pytest.warns(UserWarning, match=r"X contains .* genes that are not ."):
        exp_normalized = transformer.fit_transform(exp)
    assert isinstance(exp_normalized, pd.DataFrame)
    pd.testing.assert_frame_equal(exp_normalized, expected)

    # Validation refuses nan values
    exp.iloc[0, 0] = np.nan
    with pytest.raises(ValueError, match="Input X contains NaN."):
        transformer.fit_transform(exp)


# def test_tpm_normalization():
#     """Test TPM formula implementation.

#     TPM expressions of a minimal example have been manually computed and
#     compared against implementation.

#     Expression data:

#                         S1     S2
#     ENSG00000136807    500    400
#     ENSG00000176903  20000  19000
#     ENSG00000241490    500    600

#     Gene lengths:

#     ENSG00000136807  3000
#     ENSG00000176903  2590
#     ENSG00000241490  3101

#     Manually computed TPM (rounded):

#                          S1      S2
#     ENSG00000136807   20704   17400
#     ENSG00000176903  959266  957349
#     ENSG00000241490   20029   25250
#     """
#     genes = ["ENSG00000136807", "ENSG00000176903", "ENSG00000241490"]
#     gene_lengths = [3000, 2590, 3101]
#     expressions = [[500, 400], [20000, 19000], [500, 600]]
#     manually_computed_TPM = [
#       [20704, 17400], [959266, 957349], [20029, 25250]
#       ]

#     X = pd.DataFrame(expressions, index=genes, columns=["S1", "S2"])
#     y = pd.DataFrame(gene_lengths, index=genes, columns=["GENE_LENGTHS"])

#     TPM = tpm(X, y)
#     assert np.all(np.asarray(TPM, dtype=int) == manually_computed_TPM)

#     TPM = tpm(X.to_numpy(), y.to_numpy())
#     assert np.all(TPM.astype(int) == manually_computed_TPM)


# def test_devision_by_zero():
#     genes = ["ENSG00000136807", "ENSG00000176903", "ENSG00000241490"]
#     gene_lengths = [3000, 2590, 3101]
#     expressions = [[0], [0], [0]]
#     expected_TPM = [[0], [0], [0]]

#     X = pd.DataFrame(expressions, index=genes, columns=["S1"])
#     y = pd.DataFrame(gene_lengths, index=genes, columns=["GENE_LENGTHS"])

#     TPM = tpm(X.to_numpy(), y.to_numpy())

#     assert np.all(TPM.astype(int) == expected_TPM)


# def test_geneset_mismatch():
#     genes = ["ENSG00000136807", "ENSG00000176903", "ENSG00000241490"]
#     gene_lengths = [3000, 2590]
#     expressions = [[0], [0], [0]]
#     expected_TPM = [[0], [0]]

#     X = pd.DataFrame(expressions, index=genes, columns=["S1"])
#     y = pd.DataFrame(
#       gene_lengths, index=genes[:-1], columns=["GENE_LENGTHS"]
#       )

#     with pytest.warns(RuntimeWarning, match="Geneset mismatch"):
#         TPM = tpm(X, y)

#     assert np.all(np.asarray(TPM, dtype=int) == expected_TPM)
