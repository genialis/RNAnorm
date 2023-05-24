"""Command line interface."""
import functools
import io
import sys
from pathlib import Path
from typing import Any, Callable, Type, Union

import click
import pandas as pd

from rnanorm import CPM, CTF, CUF, FPKM, TMM, TPM, UQ

method_type = Type[Union[CPM, FPKM, TMM, TPM, UQ]]
out_type = Union[io.TextIOWrapper, Path]


class CLWrapper:
    """Input / output wrapper around normalization methods."""

    def __init__(self, method: method_type, **kwargs: Any) -> None:
        """Initialize class."""
        self.method = method
        self.kwargs = kwargs

    def handle(self, exp: pd.DataFrame, out: out_type, force: bool) -> None:
        """Parse, tranasform and write out the results.."""
        X = self.parse_exp(exp)
        normalization = self.method(**self.kwargs).set_output(transform="pandas")
        result = normalization.fit_transform(X)
        self.write_out(out, result, force)

    def parse_exp(self, exp: pd.DataFrame) -> pd.DataFrame:
        """Parse expression matrix into a DataFrame."""
        if exp.name != "<stdin>":
            if not Path(exp.name).exists():
                raise ValueError(f"File {exp.name} does not exist.")
            if not Path(exp.name).is_file():
                raise ValueError(f"File {exp.name} is not a file.")

        return pd.read_csv(exp, index_col=0)

    def write_out(self, out: out_type, result: pd.DataFrame, force: bool = False) -> None:
        """Write the resulting DataFrame to a file or stdout."""
        if out.name != "<stdout>":
            if Path(out.name).is_dir():
                raise ValueError(f"File {out.name} is a directory.")
            elif Path(out.name).exists() and not force:
                raise ValueError(
                    f"File {out.name} already exists. Use --force parameter to overwrite"
                )
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
        required=True,
        help="Compute gene lengths from this GTF file",
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


@click.command()
@common_params
def cpm(exp: pd.DataFrame, out: out_type, force: bool) -> None:
    """Compute CPM."""
    CLWrapper(CPM).handle(exp, out, force)


@click.command()
@common_params
@gtf_param
def fpkm(exp: pd.DataFrame, out: out_type, force: bool, gtf: io.TextIOWrapper) -> None:
    """Compute FPKM."""
    CLWrapper(FPKM, gtf=gtf).handle(exp, out, force)


@click.command()
@common_params
@gtf_param
def tpm(exp: pd.DataFrame, out: out_type, force: bool, gtf: io.TextIOWrapper) -> None:
    """Compute TPM."""
    CLWrapper(TPM, gtf=gtf).handle(exp, out, force)


@click.command()
@common_params
def uq(exp: pd.DataFrame, out: out_type, force: bool) -> None:
    """Compute UQ."""
    CLWrapper(UQ).handle(exp, out, force)


@click.command()
@common_params
def cuf(exp: pd.DataFrame, out: out_type, force: bool) -> None:
    """Compute CUF."""
    CLWrapper(CUF).handle(exp, out, force)


@click.command()
@common_params
@tmm_params
def tmm(exp: pd.DataFrame, out: out_type, force: bool, m_trim: float, a_trim: float) -> None:
    """Compute TMM."""
    CLWrapper(TMM, m_trim=m_trim, a_trim=a_trim).handle(exp, out, force)


@click.command()
@common_params
@tmm_params
def ctf(exp: pd.DataFrame, out: out_type, force: bool, m_trim: float, a_trim: float) -> None:
    """Compute CTF."""
    CLWrapper(CTF, m_trim=m_trim, a_trim=a_trim).handle(exp, out, force)


@click.group()
def main() -> None:
    """Define main entry point for CLI."""
    pass


# Add sub-commands to the main command
main.add_command(cpm)
main.add_command(fpkm)
main.add_command(tpm)
main.add_command(uq)
main.add_command(cuf)
main.add_command(tmm)
main.add_command(ctf)
