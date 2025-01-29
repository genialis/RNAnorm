"""Within sample normalization: CPM, FPKM and TPM."""

import warnings
from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd
from sklearn.preprocessing import FunctionTransformer
from sklearn.utils.validation import validate_data

from ..annotation import GTF
from ..typing import Numeric1D, Numeric2D
from .utils import LibrarySize


class CPM(FunctionTransformer):
    """Counts per million (CPM) normalization.

    :param allow_nan: If true, allow X to contain NaN values.

    .. rubric:: Examples

    >>> from rnanorm.datasets import load_toy_data
    >>> from rnanorm import CPM
    >>> X = load_toy_data().exp
    >>> X
              Gene_1  Gene_2  Gene_3  Gene_4  Gene_5
    Sample_1     200     300     500    2000    7000
    Sample_2     400     600    1000    4000   14000
    Sample_3     200     300     500    2000   17000
    Sample_4     200     300     500    2000    2000
    >>> CPM().set_output(transform="pandas").fit_transform(X)
               Gene_1   Gene_2    Gene_3    Gene_4    Gene_5
    Sample_1  20000.0  30000.0   50000.0  200000.0  700000.0
    Sample_2  20000.0  30000.0   50000.0  200000.0  700000.0
    Sample_3  10000.0  15000.0   25000.0  100000.0  850000.0
    Sample_4  40000.0  60000.0  100000.0  400000.0  400000.0

    """

    def __init__(self, allow_nan: bool = False) -> None:
        """Initialize class."""
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
            kwargs["ensure_all_finite"] = "allow-nan"
        # return self._validate_data(X, reset=reset, **kwargs)
        return validate_data(self, X, reset=reset, **kwargs)


class BaseNormalizationWithGTF(FunctionTransformer):
    """Normalize raw counts."""

    def __init__(
        self,
        gtf: Optional[Union[str, Path]] = None,
        gene_lengths: Optional[pd.Series] = None,
    ) -> None:
        """Initialize class."""
        if (gtf is None and gene_lengths is None) or (
            gtf is not None and gene_lengths is not None
        ):
            raise ValueError("One and only one of gtf or gene_lengths should be provided")
        self.gene_lengths = gene_lengths
        self.gtf = gtf
        super().__init__(
            func=self._func,
            validate=True,
            feature_names_out="one-to-one",
        )

    def _get_gene_lengths(self, X: Numeric2D) -> Numeric1D:
        """Check and correct X and gene lengths."""
        # gene_lengths should be a pandas.Series object with gene ID's in index
        # and gene lengths in values.
        if self.gtf:
            if not hasattr(self, "feature_names_in_"):
                raise ValueError("X should be a pandas.DataFrame object when using GTF file")
            gene_lengths = GTF(self.gtf).length
        else:
            # Validate gene lengths
            if not isinstance(self.gene_lengths, (pd.Series, np.ndarray)):
                raise ValueError("gene_lengths should be a pandas.Series or numpy.ndarray object")
            if not (self.gene_lengths > 0).all():
                raise ValueError("gene_lengths should only contain positive numbers")
            if not (self.gene_lengths - self.gene_lengths.astype(int)).sum() == 0:
                raise ValueError("gene_lengths should only contain integers")
            if isinstance(self.gene_lengths, pd.Series) and not self.gene_lengths.index.is_unique:
                raise ValueError("gene_lengths should have unique index")
            if (
                isinstance(self.gene_lengths, np.ndarray)
                and self.gene_lengths.shape[0] != X.shape[1]
            ):
                raise ValueError(
                    "gene_lengths should have same length as number of features in X."
                )

            gene_lengths = self.gene_lengths

        if hasattr(self, "feature_names_in_"):
            missing = set(self.feature_names_in_) - set(gene_lengths.index)
            if missing:
                if len(missing) == len(self.feature_names_in_):
                    container = "gtf" if self.gtf is not None else "gene_lengths"
                    raise ValueError(f"None of the genes in X are in {container}.")
                warnings.warn(
                    f"X contains {len(missing)} genes that are not in GTF "
                    f"{self.gtf}. This will result in NaN values for missing "
                    "genes in the output."
                )
            gene_lengths = gene_lengths.reindex(index=self.feature_names_in_).to_numpy()

        return gene_lengths


class FPKM(BaseNormalizationWithGTF):
    """Fragments per kilo-base million (FKPM) normalization.

    :param gtf: GTF file used when assigning counts to each gene

    .. rubric:: Examples

    >>> from rnanorm.datasets import load_toy_data
    >>> from rnanorm import FPKM
    >>> dataset = load_toy_data()
    >>> dataset.exp
              Gene_1  Gene_2  Gene_3  Gene_4  Gene_5
    Sample_1     200     300     500    2000    7000
    Sample_2     400     600    1000    4000   14000
    Sample_3     200     300     500    2000   17000
    Sample_4     200     300     500    2000    2000
    >>> fpkm = FPKM(gtf=dataset.gtf_path).set_output(transform="pandas")
    >>> fpkm.fit_transform(dataset.exp)
                Gene_1    Gene_2    Gene_3    Gene_4    Gene_5
    Sample_1  100000.0  100000.0  100000.0  200000.0  700000.0
    Sample_2  100000.0  100000.0  100000.0  200000.0  700000.0
    Sample_3   50000.0   50000.0   50000.0  100000.0  850000.0
    Sample_4  200000.0  200000.0  200000.0  400000.0  400000.0
    """

    def _func(self, X: Numeric2D) -> Numeric2D:
        gene_lengths = self._get_gene_lengths(X)
        return CPM().fit_transform(X) / gene_lengths * 1e3


class TPM(BaseNormalizationWithGTF):
    """Transcripts per kilo-base million (TPM) normalization.

    .. rubric:: Examples

    >>> from rnanorm.datasets import load_toy_data
    >>> from rnanorm import TPM
    >>> dataset = load_toy_data()
    >>> dataset.exp
              Gene_1  Gene_2  Gene_3  Gene_4  Gene_5
    Sample_1     200     300     500    2000    7000
    Sample_2     400     600    1000    4000   14000
    Sample_3     200     300     500    2000   17000
    Sample_4     200     300     500    2000    2000
    >>> tpm = TPM(gtf=dataset.gtf_path).set_output(transform="pandas")
    >>> tpm.fit_transform(dataset.exp)
                 Gene_1     Gene_2     Gene_3     Gene_4     Gene_5
    Sample_1   83333.33   83333.33   83333.33  166666.66  583333.33
    Sample_2   83333.33   83333.33   83333.33  166666.66  583333.33
    Sample_3   45454.54   45454.54   45454.54   90909.09  772727.27
    Sample_4  142857.14  142857.14  142857.14  285714.28  285714.28

    """

    def _func(self, X: Numeric2D) -> Numeric2D:
        gene_lengths = self._get_gene_lengths(X)
        # Exceptionally NaN values are allowed in case of missing genes:
        return CPM(allow_nan=True).fit_transform(X / gene_lengths * 1e3)
