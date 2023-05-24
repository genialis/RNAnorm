"""Between sample normalizations."""
import numpy as np
from sklearn.base import BaseEstimator, OneToOneFeatureMixin, TransformerMixin
from sklearn.utils.validation import check_is_fitted

from ..typing import Numeric1D, Numeric2D, Self
from .utils import LibrarySize


def remove_zero_genes(X: Numeric2D) -> Numeric2D:
    """Remove genes with zero counts in all samples.

    Also, convert to numpy.
    """
    X = np.asarray(X)
    return X[:, np.sum(X, axis=0) > 0]


def geometric_mean(x: Numeric1D) -> float:
    """Scale pd.Series array so that geometric mean of elements is 1.

    We make use of the fact that mean of the log-values is geometric mean.

    :param X: pd.Series
    :return: pd.Series scaled
    """
    return np.exp(np.mean(np.log(x)))


class UQ(OneToOneFeatureMixin, TransformerMixin, BaseEstimator):
    """Upper quartile (UQ) normalization.

    Sometimes in RNAseq a small number of genes are very highly expressed in
    some samples but not in others. This can artificially inflate library size
    and therefore (after library size normalization) cause the remaining genes
    to be considered under-sampled in those samples. Unless this effect is
    adjusted for, those genes may falsely appear to be down-regulated in that
    sample. Upper quartile is one of the approaches to correct for such
    imbalance. For more explanation on the topic check `EdgeR docs
    <https://bioconductor.org/packages/release/bioc/html/edgeR.html>`_.

    Procedure for normalization is described in `Bullard et al. 2010
    <https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-94>`_,
    but in short:

        - Use raw counts
        - Compute scaling factors
            - Remove genes that have zero counts in all samples
            - Scaling factor for each sample is the 75-th percentile
            - Rescale factors so that their geometric mean is 1
        - "Adjusted library size" = library size * normalization factors
        - Compute CPM normalization with "Adjusted library size"

    .. rubric:: Examples

    >>> from rnanorm.datasets import load_rnaseq_toy
    >>> from rnanorm import UQ
    >>> X = load_rnaseq_toy().exp
    >>> X
    array([[  200.,   300.,   500.,  2000.,  7000.],
           [  400.,   600.,  1000.,  4000., 14000.],
           [  200.,   300.,   500.,  2000., 17000.],
           [  200.,   300.,   500.,  2000.,  2000.]])
    >>> UQ().fit_transform(X)
    array([[  20000.,   30000.,   50000.,  200000.,  700000.],
           [  20000.,   30000.,   50000.,  200000.,  700000.],
           [  20000.,   30000.,   50000.,  200000., 1700000.],
           [  20000.,   30000.,   50000.,  200000.,  200000.]])

    """

    def _get_norm_factors(self, X: Numeric2D) -> Numeric1D:
        """Get UQ normalization factors (un-normalized with geometric mean).

        :param X: Expression raw count matrix (n_samples, n_features)
        """
        X = remove_zero_genes(X)
        lib_size = LibrarySize().fit_transform(X)

        # Compute upper quartile count for each sample.
        upper_quartiles = np.nanquantile(
            X,
            q=0.75,
            axis=1,
            method="nearest",
        )
        return upper_quartiles / lib_size

    def get_norm_factors(self, X: Numeric2D) -> Numeric1D:
        """Get UQ normalization factors (normalized with geometric mean).

        :param X: Expression raw count matrix (n_samples, n_features)s
        """
        check_is_fitted(self)
        factors = self._get_norm_factors(X)
        return factors / self.geometric_mean_  # type: ignore[return-value]

    def _reset(self) -> None:
        """Reset internal data-dependent state."""
        if hasattr(self, "geometric_mean_"):
            del self.geometric_mean_

    def fit(self, X: Numeric2D) -> Self:
        """Fit.

        :param X: Expression raw count matrix (n_samples, n_features)
        :return: Self
        """
        self._reset()
        X = self._validate_data(X, force_all_finite=True, reset=True)

        factors = self._get_norm_factors(X)
        self.geometric_mean_ = geometric_mean(factors)

        return self

    def transform(self, X: Numeric2D) -> Numeric2D:
        """Transform.

        :param X: Expression raw count matrix (n_samples, n_features)
        :return: Normalized expression matrix (n_samples, n_features)
        """
        check_is_fitted(self)
        X = self._validate_data(X, force_all_finite=True, reset=False)

        # Compute effective library sizes
        factors = self.get_norm_factors(X)
        effective_lib_size = LibrarySize().fit_transform(X) * factors

        # Make CPM, but with effective library size
        return X / effective_lib_size[:, np.newaxis] * 1e6


class CUF(UQ):
    """Counts adjusted with Upper quartile factors normalization.

    Procedure for normalization is described in `Johnson & Krishnan, 2022
    <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02568-9>`_,
    but in short:

        - Compute normalization factors same as in UpperQuartile
        - Divide raw counts with these factors

    .. rubric:: Examples

    >>> from rnanorm.datasets import load_rnaseq_toy
    >>> from rnanorm import CUF
    >>> X = load_rnaseq_toy().exp
    >>> X
    array([[  200.,   300.,   500.,  2000.,  7000.],
           [  400.,   600.,  1000.,  4000., 14000.],
           [  200.,   300.,   500.,  2000., 17000.],
           [  200.,   300.,   500.,  2000.,  2000.]])
    >>> CUF().fit_transform(X)
    array([[  200.,   300.,   500.,  2000.,  7000.],
           [  400.,   600.,  1000.,  4000., 14000.],
           [  400.,   600.,  1000.,  4000., 34000.],
           [  100.,   150.,   250.,  1000.,  1000.]])

    """

    def transform(self, X: Numeric2D) -> Numeric2D:
        """Transform.

        :param X: Expression raw count matrix (n_samples, n_features)
        :return: Normalized expression matrix (n_samples, n_features)
        """
        check_is_fitted(self)
        X = self._validate_data(X, force_all_finite=True, reset=False)

        # Just divide raw counts with normalization factors
        factors = self.get_norm_factors(X)
        return X / factors[:, np.newaxis]


class TMM(OneToOneFeatureMixin, TransformerMixin, BaseEstimator):
    """Trimmed mean of M-values (TMM) normalization.

    Sometime in RNAseq a small number of genes are very highly expressed in
    some samples but not in others. This can artificially inflate library size
    and therefore (after library size normalization) cause the remaining genes
    to be considered under-sampled in those samples. Unless this effect is
    adjusted for, those genes may falsely appear to be down-regulated in that
    sample. TMM is one of the approaches to correct for such imbalance. For
    more explanation on the topic check `EdgeR docs
    <https://bioconductor.org/packages/release/bioc/html/edgeR.html>`_.

    Procedure for normalization is described in `Robinson & Oshlack, 2010
    <https://doi.org/10.1186/gb-2010-11-3-r25>`_, but in short:

        - Use raw counts
        - Define the reference sample (``self.ref_``)
        - Compute scaling factors
            - Compute M values, filter by double trimming with m_trim
            - Compute A values, filter by double trimming with m_trim
            - Compute factors as weighted sum of M values
            - Factors = 2 ** factors
            - Rescale factors so that their geometric mean is 1
        - "Adjusted library size" = library size * normalization factors
        - Compute CPM normalization with "Adjusted library size"

    :param m_trim: Keep genes that are within
        (``m_trim``, 1 - ``m_trim``) percentile of M-values.
    :param a_trim: Keep genes that are within
        (``a_trim``, 1 - ``a_trim``) percentile of A-values.

    .. rubric:: Examples

    >>> from rnanorm.datasets import load_rnaseq_toy
    >>> from rnanorm import TMM
    >>> X = load_rnaseq_toy().exp
    >>> X
    array([[  200.,   300.,   500.,  2000.,  7000.],
           [  400.,   600.,  1000.,  4000., 14000.],
           [  200.,   300.,   500.,  2000., 17000.],
           [  200.,   300.,   500.,  2000.,  2000.]])
    >>> TMM().fit_transform(X)
    array([[  20000.,   30000.,   50000.,  200000.,  700000.],
           [  20000.,   30000.,   50000.,  200000.,  700000.],
           [  20000.,   30000.,   50000.,  200000., 1700000.],
           [  20000.,   30000.,   50000.,  200000.,  200000.]])

    """

    def __init__(self, m_trim: float = 0.3, a_trim: float = 0.05) -> None:
        """Initialize class."""
        self.m_trim = m_trim
        self.a_trim = a_trim

    # XXX: Could this be cached? It is called in fit and again in transform.
    # But it is hard to cache numpy arrays...
    def _get_norm_factors(self, X: Numeric2D) -> Numeric1D:
        """Get TMM normalization factors (un-normalized with geometric mean).

        :param X: Expression raw count matrix (n_samples, n_features)
        """
        X = remove_zero_genes(X)

        lib_size = LibrarySize().fit_transform(X)
        lib_size_ref = LibrarySize().fit_transform(self.ref_[np.newaxis, :])

        # Values 0 cause a lot of troubles and warnings in log / division.
        # But computing with np.nan is OK, and is handled gracefuly.
        # So convert values of 0 to np.nan early to avoid trouble.
        X[X == 0] = np.nan
        # When making 0 -> np.nan, make a copy of self.ref_, since modifying
        # self.ref_ would break library size calculation.
        ref = self.ref_.copy()
        ref[ref == 0] = np.nan

        # Equations simplify if we operate with (count / lib_size) ratio
        # Kudos for the idea to https://gitlab.com/georgy.m/conorm
        r = X / lib_size[:, np.newaxis]
        r_ref = ref / lib_size_ref

        # Gene-wise log-fold-changes, M from the paper
        m = np.log2(r) - np.log2(r_ref)
        # print(m)
        # Absolute expressions levels, A from the paper
        a = (np.log2(r) + np.log2(r_ref)) / 2
        # print(a)
        # Approximate asymptotic variances, w from the paper
        w = (1 - r) / X + (1 - r_ref) / ref
        # print(w)

        # Compute lower / upper cutoffs for M / A
        m_low, m_high = np.nanquantile(m, [self.m_trim, 1 - self.m_trim], axis=1)
        a_low, a_high = np.nanquantile(a, [self.a_trim, 1 - self.a_trim], axis=1)
        # Determine which genes to keep
        keep = np.logical_and(
            np.logical_and(m >= m_low[:, np.newaxis], m <= m_high[:, np.newaxis]),
            np.logical_and(a >= a_low[:, np.newaxis], a <= a_high[:, np.newaxis]),
        )

        # Compute weighted mean of M values for genes in keep, weights are w:
        # XXX: This seems to be a mistake in the paper: Paper says multiply
        # with w, but edgeR implementation divides with w. So do we.
        f = np.nansum(keep * m / w, axis=1) / np.nansum(keep / w, axis=1)

        # Impute nan's with 0, so that 2 ** 0 = 1
        f = np.nan_to_num(f, nan=0.0)
        factors = 2**f

        return factors

    def get_norm_factors(self, X: Numeric2D) -> Numeric1D:
        """Get UQ normalization factors (normalized with geometric mean).

        :param X: Expression raw count matrix (n_samples, n_features)s
        """
        check_is_fitted(self)
        X = self._validate_data(X, force_all_finite=True, reset=False, dtype=float)

        factors = self._get_norm_factors(X)

        return factors / self.geometric_mean_  # type: ignore[return-value]

    def _reset(self) -> None:
        """Reset internal data-dependent state."""
        if hasattr(self, "geometric_mean_"):
            del self.geometric_mean_
            del self.ref_

    def _get_ref(self, X: Numeric2D) -> Numeric1D:
        """Get reference sample."""
        f75 = UQ().fit(X).get_norm_factors(X)
        ref_index = np.argmin(abs(f75 - np.mean(f75)))
        return X[ref_index, :]

    def fit(self, X: Numeric2D) -> Self:
        """Fit.

        :param X: Expression raw count matrix (n_samples, n_features)
        :return: Self
        """
        self._reset()
        X = self._validate_data(X, force_all_finite=True, reset=True, dtype=float)
        X = remove_zero_genes(X)

        self.ref_ = self._get_ref(X)

        factors = self._get_norm_factors(X)
        self.geometric_mean_ = geometric_mean(factors)

        return self

    def transform(self, X: Numeric2D) -> Numeric2D:
        """Transform.

        :param X: Expression raw count matrix (n_samples, n_features)
        :return: Normalized expression matrix (n_samples, n_features)
        """
        # Compute effective library sizes
        factors = self.get_norm_factors(X)
        effective_lib_size = LibrarySize().fit_transform(X) * factors

        # Method ``check_is_fitted`` is not called here, since it is
        # called in self.get_norm_factors
        X = self._validate_data(X, force_all_finite=True, reset=False, dtype=float)

        # Make CPM, but with effective library size
        return X / effective_lib_size[:, np.newaxis] * 1e6


class CTF(TMM):
    """Counts adjusted with TMM factors normalization.

    Procedure for normalization is described in `Johnson & Krishnan, 2022
    <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02568-9>`_,
    but in short:

        - Compute normalization factors same as in TMM
        - Divide raw counts with these factors

    :param m_trim: Keep genes that are within
        (``m_trim``, 1 - ``m_trim``) percentile of M-values.
    :param a_trim: Keep genes that are within
        (``a_trim``, 1 - ``a_trim``) percentile of A-values.

    .. rubric:: Examples

    >>> from rnanorm.datasets import load_rnaseq_toy
    >>> from rnanorm import CTF
    >>> X = load_rnaseq_toy().exp
    >>> X
    array([[  200.,   300.,   500.,  2000.,  7000.],
           [  400.,   600.,  1000.,  4000., 14000.],
           [  200.,   300.,   500.,  2000., 17000.],
           [  200.,   300.,   500.,  2000.,  2000.]])
    >>> CTF().fit_transform(X)
    array([[  200.,   300.,   500.,  2000.,  7000.],
           [  400.,   600.,  1000.,  4000., 14000.],
           [  400.,   600.,  1000.,  4000., 34000.],
           [  100.,   150.,   250.,  1000.,  1000.]])

    """

    def transform(self, X: Numeric2D) -> Numeric2D:
        """Transform.

        :param X: Expression raw count matrix (n_samples, n_features)
        :return: Normalized expression matrix (n_samples, n_features)
        """
        # Just divide raw counts with normalization factors
        factors = self.get_norm_factors(X)

        # Method ``check_is_fitted`` is not called here, since it is
        # called in self.get_norm_factors
        X = self._validate_data(X, force_all_finite=True, reset=False, dtype=float)

        return X / factors[:, np.newaxis]
