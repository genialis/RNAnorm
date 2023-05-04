import subprocess
from io import StringIO

import pandas as pd
import pytest

from rnanorm import CPM, CTF, CUF, FPKM, TMM, TPM, UQ


@pytest.fixture
def exp():
    return pd.DataFrame(
        [
            [5, 5, 2, 5],
            [4, 1, 5, 6],
        ],
        index=pd.Series(["SAMPLE_1", "SAMPLE_2"], name="sample_id"),
        columns=["GENE_2", "GENE_1", "GENE_3", "GENE_4"],
    )


@pytest.fixture
def exp_path(tmp_path, exp):
    exp_path = tmp_path / "exp.csv"
    exp.to_csv(exp_path)
    return exp_path


@pytest.fixture
def out_path(tmp_path):
    out_path = tmp_path / "out.csv"
    yield out_path
    # Remove file after test
    out_path.unlink(missing_ok=True)


@pytest.fixture
def gtf():
    return pd.DataFrame(
        [
            [1, ".", "exon", 10, 20, ".", "+", ".", 'gene_id "GENE_1";'],
            [1, ".", "exon", 10, 20, ".", "+", ".", 'gene_id "GENE_2";'],
            [1, ".", "exon", 10, 20, ".", "+", ".", 'gene_id "GENE_3";'],
            [1, ".", "exon", 10, 20, ".", "+", ".", 'gene_id "GENE_4";'],
        ],
    )


@pytest.fixture
def gtf_path(tmp_path, gtf):
    gtf_path = tmp_path / "ann.gtf"
    gtf.to_csv(gtf_path, index=False, header=False, sep="\t")
    return gtf_path


def test_files_in_out(exp_path, out_path):
    assert not out_path.is_file()
    cp = subprocess.run(["rnanorm", "cpm", str(exp_path), "--out", str(out_path)])
    assert out_path.is_file()
    assert cp.returncode == 0


def test_stdout(exp_path, exp):
    """Output is not given, we expect to deliver via standard output."""
    cp = subprocess.run(
        ["rnanorm", "cpm", str(exp_path)],
        capture_output=True,
        text=True,
    )
    stream = StringIO(cp.stdout)
    stream.seek(0)
    pd.testing.assert_frame_equal(
        pd.read_csv(stream, index_col=0),
        CPM().set_output(transform="pandas").fit_transform(exp),
    )


def test_stdin(exp, out_path):
    """Input is not given, we expect to get input via standard input."""
    stream = StringIO()
    exp.to_csv(stream)
    stream.seek(0)

    cp = subprocess.run(
        ["rnanorm", "cpm", "--out", str(out_path)],
        input=stream.read(),
        text=True,
    )
    assert cp.returncode == 0
    pd.testing.assert_frame_equal(
        pd.read_csv(out_path, index_col=0),
        CPM().set_output(transform="pandas").fit_transform(exp),
    )


def test_fpkm(exp, exp_path, gtf_path, out_path):
    assert not out_path.is_file()

    cp = subprocess.run(
        ["rnanorm", "fpkm", str(exp_path), "--out", str(out_path), "--gtf", str(gtf_path)]
    )
    assert out_path.is_file()
    assert cp.returncode == 0
    pd.testing.assert_frame_equal(
        pd.read_csv(out_path, index_col=0),
        FPKM(gtf=gtf_path).set_output(transform="pandas").fit_transform(exp),
    )


def test_tpm(exp, exp_path, gtf_path, out_path):
    assert not out_path.is_file()

    cp = subprocess.run(
        ["rnanorm", "tpm", str(exp_path), "--out", str(out_path), "--gtf", str(gtf_path)]
    )
    assert out_path.is_file()
    assert cp.returncode == 0
    pd.testing.assert_frame_equal(
        pd.read_csv(out_path, index_col=0),
        TPM(gtf=gtf_path).set_output(transform="pandas").fit_transform(exp),
    )


def test_uq(exp, exp_path, out_path):
    assert not out_path.is_file()

    cp = subprocess.run(["rnanorm", "uq", str(exp_path), "--out", str(out_path)])
    assert out_path.is_file()
    assert cp.returncode == 0
    pd.testing.assert_frame_equal(
        pd.read_csv(out_path, index_col=0),
        UQ().set_output(transform="pandas").fit_transform(exp),
    )


def test_cuf(exp, exp_path, out_path):
    assert not out_path.is_file()

    cp = subprocess.run(["rnanorm", "cuf", str(exp_path), "--out", str(out_path)])
    assert out_path.is_file()
    assert cp.returncode == 0
    pd.testing.assert_frame_equal(
        pd.read_csv(out_path, index_col=0),
        CUF().set_output(transform="pandas").fit_transform(exp),
    )


def test_tmm(exp, exp_path, out_path):
    assert not out_path.is_file()

    cp = subprocess.run(["rnanorm", "tmm", str(exp_path), "--out", str(out_path)])
    assert out_path.is_file()
    assert cp.returncode == 0
    pd.testing.assert_frame_equal(
        pd.read_csv(out_path, index_col=0),
        TMM().set_output(transform="pandas").fit_transform(exp),
    )


def test_ctf(exp, exp_path, out_path):
    assert not out_path.is_file()

    cp = subprocess.run(["rnanorm", "ctf", str(exp_path), "--out", str(out_path)])
    assert out_path.is_file()
    assert cp.returncode == 0
    pd.testing.assert_frame_equal(
        pd.read_csv(out_path, index_col=0),
        CTF().set_output(transform="pandas").fit_transform(exp),
    )
