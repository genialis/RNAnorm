"""Between sample normalizations."""

from typing import Any, Optional

import numpy as np
import pandas as pd
from scipy.stats import gmean, rankdata, scoreatpercentile
from sklearn import config_context
from sklearn.base import BaseEstimator, OneToOneFeatureMixin, TransformerMixin
from sklearn.utils._set_output import _get_output_config
from sklearn.utils.validation import check_is_fitted, validate_data

from ..typing import Numeric1D, Numeric2D, Self
from .utils import LibrarySize


def remove_zero_genes(X: Numeric2D) -> Numeric2D:
    """Remove genes with zero count in all samples.

    Also, convert to numpy.
    """
    X = np.asarray(X)
    return X[:, np.sum(X, axis=0) > 0]


class UQ(OneToOneFeatureMixin, TransformerMixin, BaseEstimator):
    """Upper quartile (UQ) normalization.

    In an RNA-seq experiment a small fraction of genes is sometimes extremely
    overexpressed in some samples but not in others. This can artificially
    inflate library size and therefore (after library size normalization) cause
    the remaining genes to be considered under-sampled in those samples. Unless
    this effect is adjusted for, those genes may falsely appear to be
    down-regulated in that sample. Upper quartile is one of the approaches to
    correct for such imbalance. For more explanation on the topic check `EdgeR
    docs
    <https://bioconductor.org/packages/release/bioc/html/edgeR.html>`_.

    Procedure for normalization is described in `Bullard et al. 2010
    <https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-94>`_,
    but in short:

        - Use raw counts as input
        - Compute scaling factors
            - Remove genes that have zero count in all samples
            - Scaling factor is expression at the 75th percentile
            - Rescale factors so that their geometric mean is 1
        - "Adjusted library size" = library size * factor
        - Return CPM normalization with "Adjusted library size"

    This implementation is based on edgeR and has been validated to be
    identical to it to at least 10 decimal places.

    .. rubric:: Examples

    >>> from rnanorm.datasets import load_toy_data
    >>> from rnanorm import UQ
    >>> X = load_toy_data().exp
    >>> X
              Gene_1  Gene_2  Gene_3  Gene_4  Gene_5
    Sample_1     200     300     500    2000    7000
    Sample_2     400     600    1000    4000   14000
    Sample_3     200     300     500    2000   17000
    Sample_4     200     300     500    2000    2000
    >>> UQ().set_output(transform="pandas").fit_transform(X)
               Gene_1   Gene_2   Gene_3    Gene_4     Gene_5
    Sample_1  20000.0  30000.0  50000.0  200000.0   700000.0
    Sample_2  20000.0  30000.0  50000.0  200000.0   700000.0
    Sample_3  20000.0  30000.0  50000.0  200000.0  1700000.0
    Sample_4  20000.0  30000.0  50000.0  200000.0   200000.0
    """

    def _get_norm_factors(self, X: Numeric2D) -> Numeric1D:
        """Get UQ normalization factors (un-normalized with geometric mean).

        :param X: Expression raw count matrix (n_samples, n_features)
        """
        X = remove_zero_genes(X)
        # Make sure that global set_config(transform_output="pandas")
        # does not affect this method - we need numpy output here.
        lib_size = LibrarySize().set_output(transform="default").fit_transform(X)

        # Compute upper quartile count for each sample.
        # No numpy method can be used as drop-in replacement for R's quantile.
        # Scipy's method needs to be used, bit only works on 1D arrays
        upper_quartiles = np.apply_along_axis(
            func1d=scoreatpercentile,
            axis=1,
            arr=X,
            per=75,
        )
        return upper_quartiles / lib_size

    def get_norm_factors(self, X: Numeric2D) -> Numeric1D:
        """Get UQ normalization factors (normalized with geometric mean).

        :param X: Expression raw count matrix (n_samples, n_features)s
        """
        check_is_fitted(self)
        factors = self._get_norm_factors(X)
        factors = factors / self.geometric_mean_

        config = _get_output_config("transform", self)
        if config.get("dense", None) == "pandas" and isinstance(X, pd.DataFrame):
            return pd.Series(factors, index=X.index)
        return factors

    def _reset(self) -> None:
        """Reset internal data-dependent state."""
        if hasattr(self, "geometric_mean_"):
            del self.geometric_mean_

    def fit(self, X: Numeric2D, y: Optional[Numeric1D] = None, **fit_params: Any) -> Self:
        """Fit.

        :param X: Expression raw count matrix (n_samples, n_features)
        :return: Self
        """
        self._reset()
        X = validate_data(self, X, ensure_all_finite=True, reset=True)

        factors = self._get_norm_factors(X)
        self.geometric_mean_ = gmean(factors)

        return self

    def transform(self, X: Numeric2D) -> Numeric2D:
        """Transform.

        :param X: Expression raw count matrix (n_samples, n_features)
        :return: Normalized expression matrix (n_samples, n_features)
        """
        check_is_fitted(self)
        X = validate_data(self, X, ensure_all_finite=True, reset=False)

        # Compute effective library sizes
        factors = self.get_norm_factors(X)
        lib_size = LibrarySize().set_output(transform="default").fit_transform(X)
        effective_lib_size = lib_size * factors

        # Make CPM, but with effective library size
        return X / effective_lib_size[:, np.newaxis] * 1e6


class CUF(UQ):
    """Counts adjusted with Upper quartile factors normalization.

    Procedure for normalization is described in `Johnson & Krishnan, 2022
    <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02568-9>`_,
    but in short:

        - Use raw counts as input
        - Compute normalization factors the same way as in UpperQuartile
        - Divide raw counts with these factors

    .. rubric:: Examples

    >>> from rnanorm.datasets import load_toy_data
    >>> from rnanorm import CUF
    >>> X = load_toy_data().exp
    >>> X
              Gene_1  Gene_2  Gene_3  Gene_4  Gene_5
    Sample_1     200     300     500    2000    7000
    Sample_2     400     600    1000    4000   14000
    Sample_3     200     300     500    2000   17000
    Sample_4     200     300     500    2000    2000
    >>> CUF().set_output(transform="pandas").fit_transform(X)
              Gene_1  Gene_2  Gene_3  Gene_4   Gene_5
    Sample_1   200.0   300.0   500.0  2000.0   7000.0
    Sample_2   400.0   600.0  1000.0  4000.0  14000.0
    Sample_3   400.0   600.0  1000.0  4000.0  34000.0
    Sample_4   100.0   150.0   250.0  1000.0   1000.0

    """

    def transform(self, X: Numeric2D) -> Numeric2D:
        """Transform.

        :param X: Expression raw count matrix (n_samples, n_features)
        :return: Normalized expression matrix (n_samples, n_features)
        """
        check_is_fitted(self)
        X = validate_data(self, X, ensure_all_finite=True, reset=False)

        # Just divide raw counts with normalization factors
        factors = self.get_norm_factors(X)
        return X / factors[:, np.newaxis]


class TMM(OneToOneFeatureMixin, TransformerMixin, BaseEstimator):
    """Trimmed mean of M-values (TMM) normalization.

    In an RNA-seq experiment a small fraction of genes is sometimes extremely
    overexpressed in some samples but not in others . This can artificially
    inflate library size and therefore (after library size normalization) cause
    the remaining genes to be considered under-sampled in those samples. Unless
    this effect is adjusted for, those genes may falsely appear to be
    down-regulated in that sample. TMM is one of the approaches to correct for
    such imbalance. For more explanation on the topic check `EdgeR docs
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

    This implementation is based on edgeR's and is validated to be identical to
    it to at least 10 decimal places.

    :param m_trim: Keep genes that are within
        (``m_trim``, 1 - ``m_trim``) percentile of M-values.
    :param a_trim: Keep genes that are within
        (``a_trim``, 1 - ``a_trim``) percentile of A-values.

    .. rubric:: Examples

    >>> from rnanorm.datasets import load_toy_data
    >>> from rnanorm import TMM
    >>> X = load_toy_data().exp
    >>> X
              Gene_1  Gene_2  Gene_3  Gene_4  Gene_5
    Sample_1     200     300     500    2000    7000
    Sample_2     400     600    1000    4000   14000
    Sample_3     200     300     500    2000   17000
    Sample_4     200     300     500    2000    2000
    >>> TMM().set_output(transform="pandas").fit_transform(X)
               Gene_1   Gene_2   Gene_3    Gene_4     Gene_5
    Sample_1  20000.0  30000.0  50000.0  200000.0   700000.0
    Sample_2  20000.0  30000.0  50000.0  200000.0   700000.0
    Sample_3  20000.0  30000.0  50000.0  200000.0  1700000.0
    Sample_4  20000.0  30000.0  50000.0  200000.0   200000.0

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

        # ensure that output of transform will be a np.array
        with config_context(transform_output="default"):
            lib_size = LibrarySize().fit_transform(X)
            lib_size_ref = LibrarySize().fit_transform(self.ref_[np.newaxis, :])

        # Values 0 cause a lot of troubles and warnings in log / division.
        # But computing with np.nan is OK, and is handled gracefully.
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
        # Order matters: do not change from log(x/y) to log(x)-log(y)
        m = np.log2(r / r_ref)
        # Absolute expressions levels, A from the paper
        # Order matters: do not change from log(x*y) to log(x)+log(y)
        a = np.log2(r * r_ref) / 2
        # Approximate asymptotic variances, w from the paper
        w = (1 - r) / X + (1 - r_ref) / ref

        # At this point, depart from vectorization, since removal of NaN can
        # yield different lengths of rows in m / a / w in different samples It
        # would make 2D arrays with rows of diff lengths, so not numpy
        # compatible.
        f = list()
        for i in range(X.shape[0]):
            mm = m[i]
            aa = a[i]
            ww = w[i]
            # Remove NaN values
            finite = np.isfinite(mm) & np.isfinite(aa)
            mm = mm[finite]
            aa = aa[finite]
            ww = ww[finite]

            n = len(mm)
            # Compute lower / upper cutoffs for m / a values
            m_low = np.floor(n * self.m_trim) + 1
            m_high = n - m_low + 1
            a_low = np.floor(n * self.a_trim) + 1
            a_high = n - a_low + 1

            # Determine which genes to keep
            keep = np.logical_and(
                np.logical_and(rankdata(mm) >= m_low, rankdata(mm) <= m_high),
                np.logical_and(rankdata(aa) >= a_low, rankdata(aa) <= a_high),
            )

            # Compute weighted mean of M values for genes in keep, weights w:
            # XXX: This seems to be a mistake in the paper: Paper says multiply
            # with w, but edgeR implementation divides with w. So do we.
            f.append(np.nansum(keep * mm / ww) / np.nansum(keep / ww))

        f = np.asarray(f, dtype=float)  # type: ignore[assignment]
        factors = np.power(2, f)

        return factors

    def get_norm_factors(self, X: Numeric2D) -> Numeric1D:
        """Get UQ normalization factors (normalized with geometric mean).

        :param X: Expression raw count matrix (n_samples, n_features)s
        """
        check_is_fitted(self)
        X2 = validate_data(self, X, ensure_all_finite=True, reset=False, dtype=float)

        factors = self._get_norm_factors(X2)
        factors = factors / self.geometric_mean_

        config = _get_output_config("transform", self)
        if config.get("dense", None) == "pandas" and isinstance(X, pd.DataFrame):
            return pd.Series(factors, index=X.index)
        return factors

    def _reset(self) -> None:
        """Reset internal data-dependent state."""
        if hasattr(self, "geometric_mean_"):
            del self.geometric_mean_
            del self.ref_

    def _get_ref(self, X: Numeric2D) -> Numeric1D:
        """Get reference sample."""
        f75 = UQ().fit(X).get_norm_factors(X)
        ref_index = np.argmin(np.fabs(f75 - np.mean(f75)))
        return X[ref_index, :]

    def fit(self, X: Numeric2D, y: Optional[Numeric1D] = None, **fit_params: Any) -> Self:
        """Fit.

        :param X: Expression raw count matrix (n_samples, n_features)
        :return: Self
        """
        self._reset()
        X = validate_data(self, X, ensure_all_finite=True, reset=True, dtype=float)
        X = remove_zero_genes(X)

        self.ref_ = self._get_ref(X)

        factors = self._get_norm_factors(X)
        self.geometric_mean_ = gmean(factors)

        return self

    def transform(self, X: Numeric2D) -> Numeric2D:
        """Transform.

        :param X: Expression raw count matrix (n_samples, n_features)
        :return: Normalized expression matrix (n_samples, n_features)
        """
        # Compute effective library sizes
        factors = self.get_norm_factors(X)
        if isinstance(factors, pd.Series):
            factors = factors.to_numpy()

        lib_size = LibrarySize().set_output(transform="default").fit_transform(X)
        effective_lib_size = lib_size * factors

        # Method ``check_is_fitted`` is not called here, since it is
        # called in self.get_norm_factors
        X = validate_data(self, X, ensure_all_finite=True, reset=False, dtype=float)

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

    >>> from rnanorm.datasets import load_toy_data
    >>> from rnanorm import CTF
    >>> X = load_toy_data().exp
    >>> X
              Gene_1  Gene_2  Gene_3  Gene_4  Gene_5
    Sample_1     200     300     500    2000    7000
    Sample_2     400     600    1000    4000   14000
    Sample_3     200     300     500    2000   17000
    Sample_4     200     300     500    2000    2000
    >>> CTF().set_output(transform="pandas").fit_transform(X)
              Gene_1  Gene_2  Gene_3  Gene_4   Gene_5
    Sample_1   200.0   300.0   500.0  2000.0   7000.0
    Sample_2   400.0   600.0  1000.0  4000.0  14000.0
    Sample_3   400.0   600.0  1000.0  4000.0  34000.0
    Sample_4   100.0   150.0   250.0  1000.0   1000.0

    """

    def transform(self, X: Numeric2D) -> Numeric2D:
        """Transform.

        :param X: Expression raw count matrix (n_samples, n_features)
        :return: Normalized expression matrix (n_samples, n_features)
        """
        # Just divide raw counts with normalization factors
        factors = self.get_norm_factors(X)
        if isinstance(factors, pd.Series):
            factors = factors.to_numpy()

        # Method ``check_is_fitted`` is not called here, since it is
        # called in self.get_norm_factors
        X = validate_data(self, X, ensure_all_finite=True, reset=False, dtype=float)

        return X / factors[:, np.newaxis]
