"""Top-level script environment."""
import argparse
import os
import sys

import pandas as pd

from rnanorm.normalization import cpm, tpm


class ArgumentValidateException(Exception):
    """Exception for argument validation."""

    pass


class NoOutputException(Exception):
    """Exception when no output files are given."""

    pass


class InputParserException(Exception):
    """Exception when input files are improperly formatted."""

    pass


def parse_args(args):
    """Parse command line arguments."""
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
    return parser.parse_args(args)


def validate_args(args):
    """Validate that input files exist when they are needed."""
    expression_path = args.expression
    gene_lengths_path = args.gene_lengths
    tpm_output_path = args.tpm_output
    cpm_output_path = args.cpm_output

    if not os.path.exists(expression_path):
        raise ArgumentValidateException(f"Expressions file not found {expression_path}")

    if tpm_output_path:
        if gene_lengths_path is None:
            raise ArgumentValidateException("--gene-lengths must be given for --tpm-output")

        if not os.path.exists(gene_lengths_path):
            raise ArgumentValidateException(f"Gene lengths file not found {gene_lengths_path}")

    if tpm_output_path is None and cpm_output_path is None:
        raise NoOutputException


def load_expressions(expression_path):
    """Read gene expression tsv file."""
    expressions = pd.read_csv(expression_path, sep="\t")
    column_map = {col.lower(): col for col in expressions.columns}

    if "feature_id" in column_map:
        feature_id_column = column_map["feature_id"]
    else:
        raise InputParserException("FEATURE_ID column not found in expression data")

    if "gene_symbol" in column_map:
        gene_symbol_column = column_map["gene_symbol"]
        del expressions[gene_symbol_column]

    expressions.set_index(feature_id_column, inplace=True)

    return expressions


def load_gene_lengths(gene_lengths_path):
    """Read gene lengths tsv file."""
    gene_lengths = pd.read_csv(gene_lengths_path, sep="\t")
    column_map = {col.lower(): col for col in gene_lengths.columns}

    if "feature_id" in column_map:
        feature_id_column = column_map["feature_id"]
    else:
        raise InputParserException("FEATURE_ID column not found in gene_lengths data")

    gene_lengths.set_index(feature_id_column, inplace=True)

    if gene_lengths.shape[1] != 1:
        raise InputParserException("Gene lengths should have two columns: FEATURE_ID and GENE_LENGTHS")

    return gene_lengths


def main():
    """Run rnanorm command."""
    args = parse_args(sys.argv[1:])

    try:
        validate_args(args)
    except ArgumentValidateException as e:
        print(e)
        exit(1)
    except NoOutputException:
        exit(0)

    expression_path = args.expression
    gene_lengths_path = args.gene_lengths
    tpm_output_path = args.tpm_output
    cpm_output_path = args.cpm_output

    try:
        expressions = load_expressions(expression_path)

        if cpm_output_path:
            CPM = cpm(expressions)
            CPM.to_csv(cpm_output_path, sep="\t")

        if tpm_output_path:
            gene_lengths = load_gene_lengths(gene_lengths_path)
            TPM = tpm(expressions, gene_lengths)
            TPM.to_csv(tpm_output_path, sep="\t")

    except InputParserException as e:
        print(e)
        exit(1)


if __name__ == "__main__":
    main()
