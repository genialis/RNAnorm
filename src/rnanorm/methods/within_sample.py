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
    """Counts per million (CPM) normalization.

    :param allow_nan: If true, allow X to contain NaN values.

    .. rubric:: Examples

    >>> from rnanorm.datasets import load_rnaseq_toy
    >>> from rnanorm import CPM
    >>> X = load_rnaseq_toy().exp
    >>> X
    array([[  200.,   300.,   500.,  2000.,  7000.],
           [  400.,   600.,  1000.,  4000., 14000.],
           [  200.,   300.,   500.,  2000., 17000.],
           [  200.,   300.,   500.,  2000.,  2000.]])
    >>> CPM().fit_transform(X)
    array([[ 20000.,  30000.,  50000., 200000., 700000.],
           [ 20000.,  30000.,  50000., 200000., 700000.],
           [ 10000.,  15000.,  25000., 100000., 850000.],
           [ 40000.,  60000., 100000., 400000., 400000.]])

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
            kwargs["force_all_finite"] = "allow-nan"
        return self._validate_data(X, reset=reset, **kwargs)


class BaseNormalizationWithGTF(FunctionTransformer):
    """Normalize raw counts."""

    def __init__(self, gtf: Union[str, Path]) -> None:
        """Initialize class."""
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
    """Fragments per kilo-base million (FKPM) normalization.

    :param gtf: GTF file used when assigning counts to each gene

    .. rubric:: Examples

    >>> from rnanorm.datasets import load_rnaseq_toy
    >>> from rnanorm import FPKM
    >>> dataset = load_rnaseq_toy(as_frame=True)
    >>> dataset.exp
           G1     G2      G3      G4       G5
    S1  200.0  300.0   500.0  2000.0   7000.0
    S2  400.0  600.0  1000.0  4000.0  14000.0
    S3  200.0  300.0   500.0  2000.0  17000.0
    S4  200.0  300.0   500.0  2000.0   2000.0
    >>> fpkm = FPKM(gtf=dataset.gtf_path).set_output(transform="pandas")
    >>> fpkm.fit_transform(dataset.exp)
              G1        G2        G3        G4        G5
    S1  100000.0  100000.0  100000.0  200000.0  700000.0
    S2  100000.0  100000.0  100000.0  200000.0  700000.0
    S3   50000.0   50000.0   50000.0  100000.0  850000.0
    S4  200000.0  200000.0  200000.0  400000.0  400000.0

    """

    def _func(self, X: Numeric2D) -> Numeric2D:
        gene_lengths = self._get_gene_lengths(X)
        return CPM().fit_transform(X) / gene_lengths * 1e3


class TPM(BaseNormalizationWithGTF):
    """Transcripts per kilo-base million (TPM) normalization.

    .. rubric:: Examples

    >>> from rnanorm.datasets import load_rnaseq_toy
    >>> from rnanorm import TPM
    >>> dataset = load_rnaseq_toy(as_frame=True)
    >>> dataset.exp
           G1     G2      G3      G4       G5
    S1  200.0  300.0   500.0  2000.0   7000.0
    S2  400.0  600.0  1000.0  4000.0  14000.0
    S3  200.0  300.0   500.0  2000.0  17000.0
    S4  200.0  300.0   500.0  2000.0   2000.0
    >>> tpm = TPM(gtf=dataset.gtf_path).set_output(transform="pandas")
    >>> tpm.fit_transform(dataset.exp)
               G1         G2         G3         G4         G5
    S1   83333.33   83333.33   83333.33  166666.66  583333.33
    S2   83333.33   83333.33   83333.33  166666.66  583333.33
    S3   45454.54   45454.54   45454.54   90909.09  772727.27
    S4  142857.14  142857.14  142857.14  285714.28  285714.28

    """

    def _func(self, X: Numeric2D) -> Numeric2D:
        gene_lengths = self._get_gene_lengths(X)
        # Exceptionally NaN values are allowed in case of missing genes:
        return CPM(allow_nan=True).fit_transform(X / gene_lengths * 1e3)
