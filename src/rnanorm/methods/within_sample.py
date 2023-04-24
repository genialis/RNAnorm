"""Within sample normalization: CPM, FPKM and TPM."""
import warnings
from pathlib import Path
from typing import Union

import numpy as np
from sklearn.preprocessing import FunctionTransformer

from ..annotation import GTF
from ..typing import Numeric1D, Numeric2D
from .utils import LibrarySize


class CPM(FunctionTransformer):
    """Normalize raw counts with CPM (Counts per million).

    :param validate: If true, convert X to Numpy array and make
        various checks (2D matrix, non-empty, all values are finite)
        NaN values are not allowed.

    """

    def __init__(self, allow_nan: bool = False) -> None:
        self.allow_nan = allow_nan
        super().__init__(
            self._func,
            validate=True,
            feature_names_out="one-to-one",
        )

    def _func(self, X: Numeric2D) -> Numeric2D:
        """Transform.

        :param X: Expression raw count matrix (n_samples, n_features)
        :return: Normalized expression matrix (n_samples, n_features)
        """
        lib_size = LibrarySize(allow_nan=self.allow_nan).fit_transform(X)

        if hasattr(lib_size, "columns"):
            # In the case
            #   - .set_output(transform="pandas")
            #   - X is instance of pd.DataFrame
            # lib_size will be pd.DataFrame object with single
            # "Library size" column and below statement will fail with
            # InvalidIndexError. To avoid this - cast to numpy.
            lib_size = np.asarray(lib_size["Library size"])

        return X / lib_size[:, np.newaxis] * 1e6

    def _check_input(self, X: Numeric2D, *, reset: bool) -> Numeric2D:
        """Check input.

        This function is modified only so that NaN values can be
        tolerated even if self.validate=True.
        """
        kwargs = dict()
        if self.allow_nan:
            kwargs["force_all_finite"] = "allow-nan"
        return self._validate_data(X, reset=reset, **kwargs)


class BaseNormalizationWithGTF(FunctionTransformer):
    """Normalize raw counts."""

    def __init__(self, gtf: Union[str, Path]) -> None:
        self.gtf = gtf
        super().__init__(
            func=self._func,
            validate=True,
            feature_names_out="one-to-one",
        )

    def _get_gene_lengths(self, X: Numeric2D) -> Numeric1D:
        """Check and correct X and gene lengths."""
        # gene_lengths is a pandas.Series object
        gene_lengths = GTF(self.gtf).length

        missing = set(self.feature_names_in_) - set(gene_lengths.index)
        if missing:
            warnings.warn(
                f"X contains {len(missing)} genes that are not in GTF "
                f"{self.gtf}. This will result in NaN values for missing "
                "genes in the output."
            )

        return gene_lengths.reindex(index=self.feature_names_in_).to_numpy()


class FPKM(BaseNormalizationWithGTF):
    """Normalize raw counts with FKPM (Fragments per kilo-base million).

    :param build: Genome build used when assigning counts to each gene

    """

    def _func(self, X: Numeric2D) -> Numeric2D:
        gene_lengths = self._get_gene_lengths(X)
        return CPM().fit_transform(X) / gene_lengths * 1e3


class TPM(BaseNormalizationWithGTF):
    """Normalize raw counts with TPM (Transcripts per kilo-base million).

    :param build: Genome build used when assigning counts to each gene

    """

    def _func(self, X: Numeric2D) -> Numeric2D:
        gene_lengths = self._get_gene_lengths(X)
        # Exceptionally NaN values are allowed in case of missing genes:
        return CPM(allow_nan=True).fit_transform(X / gene_lengths * 1e3)
