"""TPM normalization"""
import argparse
import os
import warnings

import numpy as np
import pandas as pd


def _tpm_ndarray(X, y):
    """TPM normalization.

    :type X: Numpy ndarray
    :type y: Numpy ndarray
    """
    assert isinstance(X, np.ndarray)
    assert isinstance(y, np.ndarray)
    assert X.shape[0] == y.shape[0]
    assert y.shape[1] == 1
    assert np.min(X) >= 0.0  # Gene counts must be non-negative

    A = X / y
    sumA = A.sum(axis=0)

    with np.errstate(invalid="ignore"):  # Ignore warnings of division by 0
        TPM = A / sumA * 1e6

        # Samples with zeros for all genes get nan but should be 0.0
        np.nan_to_num(TPM, copy=False)

    return TPM


def tpm(X, y):
    """Transcript per kilobase million (TPM).

    A = readsMappedToGene / geneLength
    TPM = A / SUM(A) * 1e6

    :type X: 2-D array_like
    :type y: 1-D array_like
    """
    if isinstance(X, pd.DataFrame) and isinstance(y, pd.DataFrame):
        common_genes = X.index.intersection(y.index)
        ncommon = len(common_genes)
        if ncommon != X.shape[0] or ncommon != y.shape[0]:
            warnings.warn(
                f"Geneset mismatch between expressions ({X.shape[0]}) and gene "
                f"lengths ({y.shape[0]}). Using intersection ({ncommon})...",
                RuntimeWarning,
                stacklevel=2,
            )

        X = X.loc[common_genes]
        y = y.loc[common_genes]

        X_ = np.asarray(X, dtype=np.float64)
        y_ = np.asarray(y, dtype=np.float64)
        TPM = _tpm_ndarray(X_, y_)

        return pd.DataFrame(TPM, index=X.index, columns=X.columns)
    else:
        X = np.asarray(X, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        return _tpm_ndarray(X, y)


def cpm(X):
    """Counts per million (CPM).

    CPM = readsMappedToGene / totalNumReads * 1e6

    :type X: 2-D array_like
    """
    if isinstance(X, pd.DataFrame):
        assert X.to_numpy().min() >= 0.0  # Gene counts must be non-negative
    else:
        # Cast non-Pandas array_like objects to Numpy
        X = np.asarray(X, dtype=np.float64)
        assert np.min(X) >= 0.0  # Gene counts must be non-negative

    sumX = X.sum(axis=0)

    with np.errstate(invalid="ignore"):  # Ignore warnings of division by 0
        CPM = X / sumX * 1e6

        # Samples with zeros for all genes get nan but should be 0.0
        np.nan_to_num(CPM, copy=False)

    return CPM


def main():
    parser = argparse.ArgumentParser(
        description="""TPM normalization.

    The gene expressions file should include genes in rows and samples in columns.
    The gene ID column should be named FEATURE_ID.

    The gene lengths file should have two columns, FEATURE_ID and GENE_LENGTHS.

    Gene IDs in expressions file should match the gene IDs in gene lengths file.

    """
    )
    parser.add_argument(
        "expressions_file",
        type=str,
        help="tab-delimited file with gene expression data (genes in rows, samples in cols)",
    )
    parser.add_argument("gene_lengths_file", type=str, help="tab-delimited file with gene lengths")
    parser.add_argument("--output", "-o", type=str, help="output file name")

    args = parser.parse_args()
    expressions_path = args.expressions_file
    gene_lengths_path = args.gene_lengths_file
    output_path = args.output

    if not os.path.exists(expressions_path):
        print(f"File not found {expressions_path}")
        exit(1)

    if not os.path.exists(gene_lengths_path):
        print(f"File not found {gene_lengths_path}")
        exit(1)

    expressions = pd.read_csv(expressions_path, sep="\t")

    column_map = {col.lower(): col for col in expressions.columns}
    if "feature_id" in column_map:
        feature_id_column = column_map["feature_id"]
    else:
        print("FEATURE_ID column not found in expression data")
        exit(1)

    if "gene_symbol" in column_map:
        gene_symbol_column = column_map["gene_symbol"]
        del expressions[gene_symbol_column]

    expressions.set_index(feature_id_column, inplace=True)

    gene_lengths = pd.read_csv(gene_lengths_path, sep="\t")
    column_map = {col.lower(): col for col in gene_lengths.columns}
    if "feature_id" in column_map:
        feature_id_column = column_map["feature_id"]
    else:
        print("FEATURE_ID column not found in gene_lengths data")
        exit(1)

    gene_lengths.set_index(feature_id_column, inplace=True)

    if gene_lengths.shape[1] != 1:
        print("Gene lengths should have two columns: FEATURE_ID and GENE_LENGTHS")
        exit(1)

    TPM = tpm(expressions, gene_lengths)

    if output_path is None:
        print(TPM.to_csv(sep="\t"))
    else:
        TPM.to_csv(output_path, sep="\t")
