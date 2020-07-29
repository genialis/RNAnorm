"""Command-line commands."""
import argparse
import os

import pandas as pd

from .normalization import cpm, tpm


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
        "expression", type=str, help="tab-delimited file with gene expression data (genes in rows, samples in cols)",
    )
    parser.add_argument("--gene-lengths", type=str, help="tab-delimited file with gene lengths")
    parser.add_argument("--tpm-output", type=str, help="TPM output file name")
    parser.add_argument("--cpm-output", type=str, help="CPM output file name")

    args = parser.parse_args()
    expression_path = args.expression
    gene_lengths_path = args.gene_lengths
    tpm_output_path = args.tpm_output
    cpm_output_path = args.cpm_output

    if not os.path.exists(expression_path):
        print(f"Expressions file not found {expression_path}")
        exit(1)

    if tpm_output_path:
        if gene_lengths_path is None:
            print("--gene-lengths must be given for --tpm-output")
            exit(1)

        if not os.path.exists(gene_lengths_path):
            print(f"Gene lengths file not found {gene_lengths_path}")
            exit(1)

    if tpm_output_path is None and cpm_output_path is None:
        exit()

    expressions = pd.read_csv(expression_path, sep="\t")

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

    if cpm_output_path:
        CPM = cpm(expressions)
        CPM.to_csv(cpm_output_path, sep="\t")

    if tpm_output_path:
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
        TPM.to_csv(tpm_output_path, sep="\t")
