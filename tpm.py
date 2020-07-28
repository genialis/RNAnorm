"""TPM normalization"""
import argparse
import os

import numpy as np
import pandas as pd


__version__ = "1.0.0"

parser = argparse.ArgumentParser(description='TPM normalization')
parser.add_argument('expressions_file', type=str, help='tab file with gene expression data (genes in rows, samples in cols)')
parser.add_argument('gene_lengths_file', type=str, help='tab file with gene lengths')


def tpm_normalization(X, l):
    assert X.shape[0] == l.shape[0]
    assert l.shape[1] == 1
    assert isinstance(X, np.ndarray)
    assert isinstance(l, np.ndarray)

    A = X / l
    sumA = A.sum(axis=1)
    sumA = sumA[:,None]  # Create column from vector
    TPM = 1e6 * A / sumA
    return TPM


def format_inputs(X, l):
    assert isinstance(X, pd.DataFrame)
    assert isinstance(l, pd.DataFrame)
    X = X.to_numpy()
    l = l.to_numpy()
    return X, l


def main():
    args = parser.parse_args()
    expressions_path = args.expressions_file
    gene_lengths_path = args.gene_lengths_file

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
        print("FEATURE_ID column not found in expression data")
        exit(1)

    gene_lengths.set_index(feature_id_column, inplace = True)

    X, l = format_inputs(expressions, gene_lengths)
    tpm_normalization(X, l)


if __name__ == "__name__":
    main()
