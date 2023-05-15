import subprocess
from io import StringIO

import pandas as pd
import pytest

from rnanorm import CPM, CTF, CUF, FPKM, TMM, TPM, UQ


@pytest.fixture
def exp_path(tmp_path, exp):
    exp_path = tmp_path / "exp.csv"

    # This si to prevent warnings for missing gene in TPM / FPKM
    exp.drop(columns="G5", inplace=True)

    exp.to_csv(exp_path)
    return exp_path


@pytest.fixture
def out_path(tmp_path):
    out_path = tmp_path / "out.csv"
    yield out_path
    # Remove file after test
    out_path.unlink(missing_ok=True)


def test_files_in_out(exp_path, out_path):
    assert not out_path.is_file()
    cp = subprocess.run(["rnanorm", "cpm", str(exp_path), "--out", str(out_path)])
    assert out_path.is_file()
    assert cp.returncode == 0


def test_stdout(exp_path, exp):
    """Parameter --out is not given, deliver via standard output."""
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
    """Input is not given, get input via standard input."""
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


def test_fpkm(exp, exp_path, gtf_file, out_path):
    assert not out_path.is_file()

    cp = subprocess.run(
        ["rnanorm", "fpkm", str(exp_path), "--out", str(out_path), "--gtf", str(gtf_file)]
    )
    assert out_path.is_file()
    assert cp.returncode == 0
    pd.testing.assert_frame_equal(
        pd.read_csv(out_path, index_col=0),
        FPKM(gtf=gtf_file).set_output(transform="pandas").fit_transform(exp),
    )


def test_tpm(exp, exp_path, gtf_file, out_path):
    assert not out_path.is_file()

    cp = subprocess.run(
        ["rnanorm", "tpm", str(exp_path), "--out", str(out_path), "--gtf", str(gtf_file)]
    )
    assert out_path.is_file()
    assert cp.returncode == 0
    pd.testing.assert_frame_equal(
        pd.read_csv(out_path, index_col=0),
        TPM(gtf=gtf_file).set_output(transform="pandas").fit_transform(exp),
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
