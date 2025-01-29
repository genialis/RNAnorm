"""Command line interface."""

import functools
import io
import sys
from pathlib import Path
from typing import Any, Callable, Optional, Type, Union

import click
import pandas as pd

from rnanorm import CPM, CTF, CUF, FPKM, TMM, TPM, UQ

method_type = Type[Union[CPM, FPKM, TMM, TPM, UQ]]
file_type = Union[io.TextIOWrapper, Path]


class CLWrapper:
    """Input / output wrapper around normalization methods."""

    def __init__(self, method: method_type, **kwargs: Any) -> None:
        """Initialize class."""
        self.method = method
        self.kwargs = kwargs

    def handle(self, exp: file_type, out: file_type, force: bool) -> None:
        """Parse, transform and write out the results."""
        X = self.parse_exp(exp)

        # Turn ``--gene-lengths`` (CLI option, a file object) into
        # ``gene_lengths`` (TPM / FPKM method argument, a Series)
        if self.kwargs.get("gene_lengths", None) is not None:
            self.kwargs["gene_lengths"] = self.parse_gene_lengths(self.kwargs["gene_lengths"])

        normalization = self.method(**self.kwargs).set_output(transform="pandas")
        result = normalization.fit_transform(X)
        self.write_out(out, result, force)

    def parse_exp(self, exp: file_type) -> pd.DataFrame:
        """Parse expression matrix into a DataFrame."""
        if exp.name != "<stdin>":
            if not Path(exp.name).exists():
                raise ValueError(f"File {exp.name} does not exist.")
            if not Path(exp.name).is_file():
                raise ValueError(f"File {exp.name} is not a file.")

        return pd.read_csv(exp, index_col=0)

    def parse_gene_lengths(self, gene_lengths_file: file_type) -> pd.Series:
        """Parse gene lengths file into a Series object."""
        if not Path(gene_lengths_file.name).exists():
            raise ValueError(f"File {gene_lengths_file.name} does not exist.")
        if not Path(gene_lengths_file.name).is_file():
            raise ValueError(f"File {gene_lengths_file.name} is not a file.")

        return pd.read_csv(gene_lengths_file, index_col=0).iloc[:, 0]

    def write_out(self, out: file_type, result: pd.DataFrame, force: bool = False) -> None:
        """Write the resulting DataFrame to a file or stdout."""
        if out.name != "<stdout>":
            if Path(out.name).is_dir():
                raise ValueError(f"File {out.name} is a directory.")
            elif Path(out.name).exists() and not force:
                raise ValueError(f"File {out.name} already exists. Use --force flag to overwrite.")
            elif not Path(out.name).parent.exists():
                Path(out.name).parent.mkdir(parents=True)

        # result.index.name = "sample_id"
        result.to_csv(out)


def common_params(func: Callable[..., Any]) -> Callable[..., Any]:
    """Set common parameters for all normalization methods."""

    @click.argument("exp", type=click.File("r"), default=sys.stdin)
    @click.option(
        "--out",
        type=click.File("w"),
        default=sys.stdout,
        help="Output results in this file instead of stdout",
    )
    @click.option("--force", is_flag=True, help="Overwrite already existing output file")
    @functools.wraps(func)
    def wrapper(*args: Any, **kwargs: Any) -> Callable[..., Any]:
        return func(*args, **kwargs)

    return wrapper


def gtf_param(func: Callable[..., Any]) -> Callable[..., Any]:
    """Set parameters for normalization methods that require a GTF file."""

    @click.option(
        "--gtf",
        type=click.File("r"),
        required=False,
        help="Compute gene lengths from this GTF file",
    )
    @click.option(
        "--gene-lengths",
        type=click.File("r"),
        required=False,
        help="File with gene lengths",
    )
    @functools.wraps(func)
    def wrapper(*args: Any, **kwargs: Any) -> Callable[..., Any]:
        return func(*args, **kwargs)

    return wrapper


def tmm_params(func: Callable[..., Any]) -> Callable[..., Any]:
    """Set parameters for TMM and CTF normalization."""

    @click.option("--m_trim", default=0.3, help="Two sided cutoff for M-values")
    @click.option("--a_trim", default=0.05, help="Two sided cutoff for A-values")
    @functools.wraps(func)
    def wrapper(*args: Any, **kwargs: Any) -> Callable[..., Any]:
        return func(*args, **kwargs)

    return wrapper


@click.command(short_help="Counts per million")
@common_params
def cpm(exp: pd.DataFrame, out: file_type, force: bool) -> None:
    """Compute CPM."""
    CLWrapper(CPM).handle(exp, out, force)


@click.command(short_help="Fragments per kilo-base million")
@common_params
@gtf_param
def fpkm(
    exp: pd.DataFrame,
    out: file_type,
    force: bool,
    gtf: Optional[io.TextIOWrapper] = None,
    gene_lengths: Optional[io.TextIOWrapper] = None,
) -> None:
    """Compute FPKM."""
    CLWrapper(FPKM, gtf=gtf, gene_lengths=gene_lengths).handle(exp, out, force)


@click.command(short_help="Transcripts per million")
@common_params
@gtf_param
def tpm(
    exp: pd.DataFrame,
    out: file_type,
    force: bool,
    gtf: Optional[io.TextIOWrapper] = None,
    gene_lengths: Optional[io.TextIOWrapper] = None,
) -> None:
    """Compute TPM."""
    CLWrapper(TPM, gtf=gtf, gene_lengths=gene_lengths).handle(exp, out, force)


@click.command(short_help="Upper quartile")
@common_params
def uq(exp: pd.DataFrame, out: file_type, force: bool) -> None:
    """Compute UQ."""
    CLWrapper(UQ).handle(exp, out, force)


@click.command(short_help="Counts adjusted with UQ factors")
@common_params
def cuf(exp: pd.DataFrame, out: file_type, force: bool) -> None:
    """Compute CUF."""
    CLWrapper(CUF).handle(exp, out, force)


@click.command(short_help="Trimmed mean of M-values")
@common_params
@tmm_params
def tmm(exp: pd.DataFrame, out: file_type, force: bool, m_trim: float, a_trim: float) -> None:
    """Compute TMM."""
    CLWrapper(TMM, m_trim=m_trim, a_trim=a_trim).handle(exp, out, force)


@click.command(short_help="Counts adjusted with TMM factors")
@common_params
@tmm_params
def ctf(exp: pd.DataFrame, out: file_type, force: bool, m_trim: float, a_trim: float) -> None:
    """Compute CTF."""
    CLWrapper(CTF, m_trim=m_trim, a_trim=a_trim).handle(exp, out, force)


@click.group()
def main() -> None:
    """Common RNA-seq normalization methods."""
    pass


# Add sub-commands to the main command
main.add_command(cpm)
main.add_command(fpkm)
main.add_command(tpm)
main.add_command(uq)
main.add_command(cuf)
main.add_command(tmm)
main.add_command(ctf)
