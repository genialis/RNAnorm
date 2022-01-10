import csv

import numpy as np
import pandas as pd
import pytest

from rnanorm.__main__ import (
    ArgumentValidateException,
    InputParserException,
    load_expressions,
    load_gene_lengths,
    parse_args,
    validate_args,
)
from rnanorm.normalization import cpm, fpkm, tpm


@pytest.fixture(scope="session")
def create_expression_files(tmpdir_factory):
    data_dir = tmpdir_factory.mktemp("data")

    genes = ["ENSG00000136807", "ENSG00000176903", "ENSG00000241490"]
    gene_lengths = [3000, 2590, 3101]
    expressions = [[500, 400], [20000, 19000], [500, 600]]

    expression_path = data_dir.join("expr.tsv")
    with open(expression_path, "w") as fp:
        writer = csv.writer(fp, delimiter="\t")
        writer.writerow(["FEATURE_ID", "S1", "S2"])
        for gene, expr in zip(genes, expressions):
            writer.writerow([gene] + expr)

    gene_lengths_path = data_dir.join("lengths.tsv")
    with open(gene_lengths_path, "w") as fp:
        writer = csv.writer(fp, delimiter="\t")
        writer.writerow(["FEATURE_ID", "GENE_LENGTHS"])
        for gene, length in zip(genes, gene_lengths):
            writer.writerow([gene, length])

    return expression_path, gene_lengths_path


@pytest.fixture(scope="session")
def create_output_files(tmpdir_factory):
    output_dir = tmpdir_factory.mktemp("output")
    tpm_output_path = output_dir.join("expr.tpm.tsv")
    fpkm_output_path = output_dir.join("expr.fpkm.tsv")
    cpm_output_path = output_dir.join("expr.cpm.tsv")

    return tpm_output_path, fpkm_output_path, cpm_output_path


def test_parse_args():
    """Test argument parser."""
    cmdline = [
        "--gene-lengths=lengths.tsv",
        "--cpm-output=expr.cpm.tsv",
        "--tpm-output=expr.tpm.tsv",
        "--fpkm-output=expr.fpkm.tsv",
        "expr.tsv",
    ]
    args = parse_args(cmdline)

    checks = [
        args.expression == "expr.tsv",
        args.gene_lengths == "lengths.tsv",
        args.tpm_output == "expr.tpm.tsv",
        args.fpkm_output == "expr.fpkm.tsv",
        args.cpm_output == "expr.cpm.tsv",
    ]
    assert all(checks)


def test_validate_args():
    """Test for existence check for gene expression file.

    Test if validation function raises an exception when gene expression file is missing.
    """
    with pytest.raises(ArgumentValidateException):
        args = parse_args(["mock_filename"])
        validate_args(args)


def test_validate_tpm_args(create_expression_files, create_output_files):
    """Test for existence of gene lengths file when TPM normalization type is used."""
    expression_path = create_expression_files[0]
    tpm_output_path = create_output_files[0]

    cmdline = [f"--tpm-output={tpm_output_path}", f"{expression_path}"]

    with pytest.raises(ArgumentValidateException):
        args = parse_args(cmdline)
        validate_args(args)


def test_validate_fpkm_args(create_expression_files, create_output_files):
    """Test for existence of gene lengths file when FPKM normalization type is used."""
    expression_path = create_expression_files[0]
    fpkm_output_path = create_output_files[1]

    cmdline = [f"--fpkm-output={fpkm_output_path}", f"{expression_path}"]

    with pytest.raises(ArgumentValidateException):
        args = parse_args(cmdline)
        validate_args(args)


def test_normalization_output_success(create_expression_files, create_output_files):
    """Test TPM, FPKM and CPM output.

    This test implements the contents of main function and checks if the
    normalization outputs are correctly written to files.
    """
    expression_path, gene_lengths_path = create_expression_files
    tpm_output_path, fpkm_output_path, cpm_output_path = create_output_files

    cmdline = [
        f"--gene-lengths={gene_lengths_path}",
        f"--cpm-output={cpm_output_path}",
        f"--tpm-output={tpm_output_path}",
        f"--fpkm-output={fpkm_output_path}",
        f"{expression_path}",
    ]

    args = parse_args(cmdline)

    # Test if validation passes when files exist
    validate_args(args)

    # Test if tpm and cpm can be computed and written to files
    expressions = load_expressions(expression_path)

    if cpm_output_path:
        CPM = cpm(expressions)
        CPM.to_csv(cpm_output_path, sep="\t")

    if tpm_output_path:
        gene_lengths = load_gene_lengths(gene_lengths_path)
        TPM = tpm(expressions, gene_lengths)
        TPM.to_csv(tpm_output_path, sep="\t")
    if fpkm_output_path:
        gene_lengths = load_gene_lengths(gene_lengths_path)
        FPKM = fpkm(expressions, gene_lengths)
        FPKM.to_csv(fpkm_output_path, sep="\t")

    # Re-read the values in files and check if they are correct
    manually_computed_CPM = [[23809, 20000], [952380, 950000], [23809, 30000]]
    CPM = pd.DataFrame(pd.read_csv(cpm_output_path, sep="\t"), columns=["S1", "S2"])
    assert np.all(np.asarray(CPM, dtype=int) == manually_computed_CPM)

    manually_computed_TPM = [[20704, 17400], [959266, 957349], [20029, 25250]]
    TPM = pd.DataFrame(pd.read_csv(tpm_output_path, sep="\t"), columns=["S1", "S2"])
    assert np.all(np.asarray(TPM, dtype=int) == manually_computed_TPM)

    manually_computed_FPKM = [[7936, 6666], [367714, 366795], [7678, 9674]]
    FPKM = pd.DataFrame(pd.read_csv(fpkm_output_path, sep="\t"), columns=["S1", "S2"])
    assert np.all(np.asarray(FPKM, dtype=int) == manually_computed_FPKM)


def test_header_validation(create_expression_files, create_output_files):
    """Test header validation in expression file.

    Test if expression load function check for existence of FEATURE_ID column.
    """
    expression_path = create_expression_files[0]
    expression_path_headerless = expression_path + ".nohead"

    # Remove header, use different file to prevent changing session-wide expressions file
    with open(expression_path, "r") as fp:
        data = fp.read().splitlines(True)
    with open(expression_path_headerless, "w") as fp:
        for line in data:
            if "FEATURE_ID" in line:
                line = line.replace("FEATURE_ID", "MOCK_ID")
            fp.write(line)

    with pytest.raises(InputParserException):
        load_expressions(expression_path_headerless)
