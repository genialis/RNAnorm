"""Top-level script environment."""
import argparse
import os
import sys

import pandas as pd

from rnanorm.annotation import union_exon_lengths
from rnanorm.normalization import cpm, fpkm, quantile, tpm


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
        "expression",
        type=str,
        help="tab-delimited file with gene expression data (genes in rows, samples in cols)",
    )
    parser.add_argument("--gene-lengths", type=str, help="tab-delimited file with gene lengths")
    parser.add_argument("--annotation", type=str, help="Annotation file in GTF format")
    parser.add_argument("--gene-id-attr", type=str, help="Gene ID attribute for annotation file", default="gene_id")
    parser.add_argument("--tpm-output", type=str, help="TPM output file name")
    parser.add_argument("--fpkm-output", type=str, help="FPKM output file name")
    parser.add_argument("--cpm-output", type=str, help="CPM output file name")
    parser.add_argument("--quantile-output", type=str, help="Quantile-normalized expression output file name")
    return parser.parse_args(args)


def validate_args(args):
    """Validate that input files exist when they are needed."""
    expression_path = args.expression
    gene_lengths_path = args.gene_lengths
    annotation_path = args.annotation
    tpm_output_path = args.tpm_output
    fpkm_output_path = args.fpkm_output
    cpm_output_path = args.cpm_output
    quantile_output_path = args.quantile_output

    if not os.path.exists(expression_path):
        raise ArgumentValidateException(f"Expressions file not found {expression_path}")

    if tpm_output_path or fpkm_output_path:
        # 1. Validate argument consistency
        if gene_lengths_path is None and annotation_path is None:
            raise ArgumentValidateException("--gene-lengths or --annotation must be given")

        if gene_lengths_path is not None and annotation_path is not None:
            raise ArgumentValidateException("Only one of --gene-lengths or --annotation must be given")

        # 2. Validate that files actually exist.
        if gene_lengths_path and not os.path.exists(gene_lengths_path):
            raise ArgumentValidateException(f"Gene lengths file not found {gene_lengths_path}")
        if annotation_path and not os.path.exists(annotation_path):
            raise ArgumentValidateException(f"Annotation file not found {annotation_path}")

    if all([path is None for path in [tpm_output_path, fpkm_output_path, cpm_output_path, quantile_output_path]]):
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
    annotation_path = args.annotation
    gene_id_attr = args.gene_id_attr
    tpm_output_path = args.tpm_output
    fpkm_output_path = args.fpkm_output
    cpm_output_path = args.cpm_output
    quantile_output_path = args.quantile_output

    try:
        expressions = load_expressions(expression_path)

        if cpm_output_path:
            CPM = cpm(expressions)
            CPM.to_csv(cpm_output_path, sep="\t")

        if tpm_output_path or fpkm_output_path:
            if gene_lengths_path:
                gene_lengths = load_gene_lengths(gene_lengths_path)
            elif annotation_path:
                gene_lengths = union_exon_lengths(annotation_path, gene_id_attr)

            if tpm_output_path:
                TPM = tpm(expressions, gene_lengths)
                TPM.to_csv(tpm_output_path, sep="\t")

            if fpkm_output_path:
                FPKM = fpkm(expressions, gene_lengths)
                FPKM.to_csv(fpkm_output_path, sep="\t")

        if quantile_output_path:
            QT = quantile(expressions)
            QT.to_csv(quantile_output_path, sep="\t")

    except InputParserException as e:
        print(e)
        exit(1)


if __name__ == "__main__":
    main()
