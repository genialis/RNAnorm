import pandas as pd

from rnanorm.normalization import tpm_normalization


def test_tpm_normalization():
    expressions = [[500, 400], [20000, 19000], [500, 600]]
    genes = ["ENSG00000136807", "ENSG00000176903", "ENSG00000241490"]
    gene_lengths = [3000, 2590, 3101]
    X = pd.DataFrame(expressions, index=genes, columns=["S1", "S2"])
    y = pd.DataFrame(gene_lengths, index=genes, columns=["GENE_LENGTHS"])
    s11 = X.loc["ENSG00000136807"]["S1"]
    s12 = X.loc["ENSG00000176903"]["S1"]
    s13 = X.loc["ENSG00000241490"]["S1"]
    l1 = y.loc["ENSG00000136807"]["GENE_LENGTHS"]
    l2 = y.loc["ENSG00000176903"]["GENE_LENGTHS"]
    l3 = y.loc["ENSG00000241490"]["GENE_LENGTHS"]
    ##################################################
    # New implementation according with one web source
    a11 = s11 / l1
    a12 = s12 / l2
    a13 = s13 / l3
    a = a11 + a12 + a13
    g11 = a11 / a * 1e6
    g12 = a12 / a * 1e6
    g13 = a13 / a * 1e6
    print(g11)
    print(g12)
    print(g13)
    print()
    ################################################
    # Janez'z formula picked from another web source
    rpk1 = s11 * 1e3 / l1
    rpk2 = s12 * 1e3 / l2
    rpk3 = s13 * 1e3 / l3
    rpk = rpk1 + rpk2 + rpk3
    g11 = rpk1 / rpk * 1e6
    g12 = rpk2 / rpk * 1e6
    g13 = rpk3 / rpk * 1e6
    print(g11)
    print(g12)
    print(g13)
    print()
    TPM = tpm_normalization(X.to_numpy(), y.to_numpy())
    print(TPM)
