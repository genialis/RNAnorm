"""RNA-seq normalization methods."""
import warnings

import numpy as np
import pandas as pd
from scipy.stats import rankdata


def _tpm_ndarray(X, y):
    """Normalize Numpy ndarray expression counts to TPM.

    :type X: Numpy ndarray
    :type y: Numpy ndarray
    """
    assert isinstance(X, np.ndarray)
    assert isinstance(y, np.ndarray)
    assert X.shape[0] == y.shape[0]
    assert y.shape[1] == 1
    assert np.min(X) >= 0.0  # Gene counts must be non-negative

    A = X / y
    sumA = A.sum(axis=0)

    with np.errstate(invalid="ignore"):  # Ignore warnings of division by 0
        TPM = 1e6 * A / sumA

        # Samples with zeros for all genes get nan but should be 0.0
        np.nan_to_num(TPM, copy=False)

    return TPM


def tpm(X, y):
    """Normalize expression counts to Transcript per kilobase million (TPM).

    A = readsMappedToGene / geneLength
    TPM = A / SUM(A) * 1e6

    :type X: 2-D array_like
    :type y: 1-D array_like
    """
    if isinstance(X, pd.DataFrame) and isinstance(y, pd.DataFrame):
        common_genes = X.index.intersection(y.index)
        ncommon = len(common_genes)
        if ncommon != X.shape[0] or ncommon != y.shape[0]:
            warnings.warn(
                f"Geneset mismatch between expressions ({X.shape[0]}) and gene "
                f"lengths ({y.shape[0]}). Using intersection ({ncommon})...",
                RuntimeWarning,
                stacklevel=2,
            )

        X = X.loc[common_genes]
        y = y.loc[common_genes]

        X_ = np.asarray(X, dtype=np.float64)
        y_ = np.asarray(y, dtype=np.float64)
        TPM = _tpm_ndarray(X_, y_)

        return pd.DataFrame(TPM, index=X.index, columns=X.columns)
    else:
        X = np.asarray(X, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        return _tpm_ndarray(X, y)


def _fpkm_ndarray(X, y):
    """Normalize Numpy ndarray expression counts to FPKM.

    :type X: Numpy ndarray
    :type y: Numpy ndarray
    """
    assert isinstance(X, np.ndarray)
    assert isinstance(y, np.ndarray)
    assert X.shape[0] == y.shape[0]
    assert y.shape[1] == 1
    assert np.min(X) >= 0.0  # Gene counts must be non-negative

    total_sample_reads = X.sum(axis=0) / 10**6

    with np.errstate(invalid="ignore"):  # Ignore warnings of division by 0
        rpm = X / total_sample_reads
        # Samples with zeros for all genes get nan but should be 0.0
        np.nan_to_num(rpm, copy=False)

    fpkm = rpm / y * 1000

    return fpkm


def fpkm(X, y):
    """Normalize expression counts to Fragments per kilobase million (FPKM).

    RPM = expressionCount / sumReadsInSample * 1e6
    FPKM = RPM / geneLengthKb

    :type X: 2-D array_like
    :type y: 1-D array_like
    """
    if isinstance(X, pd.DataFrame) and isinstance(y, pd.DataFrame):
        common_genes = X.index.intersection(y.index)
        ncommon = len(common_genes)
        if ncommon != X.shape[0] or ncommon != y.shape[0]:
            warnings.warn(
                f"Geneset mismatch between expressions ({X.shape[0]}) and gene "
                f"lengths ({y.shape[0]}). Using intersection ({ncommon})...",
                RuntimeWarning,
                stacklevel=2,
            )

        X = X.loc[common_genes]
        y = y.loc[common_genes]

        X_ = np.asarray(X, dtype=np.float64)
        y_ = np.asarray(y, dtype=np.float64)
        FPKM = _fpkm_ndarray(X_, y_)

        return pd.DataFrame(FPKM, index=X.index, columns=X.columns)
    else:
        X = np.asarray(X, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        return _fpkm_ndarray(X, y)


def cpm(X):
    """Normalize expression counts to Counts per million (CPM).

    CPM = readsMappedToGene / totalNumReads * 1e6

    :type X: 2-D array_like
    """
    if isinstance(X, pd.DataFrame):
        assert X.to_numpy().min() >= 0.0  # Gene counts must be non-negative
    else:
        # Cast non-Pandas array_like objects to Numpy
        X = np.asarray(X, dtype=np.float64)
        assert np.min(X) >= 0.0  # Gene counts must be non-negative

    sumX = X.sum(axis=0)

    with np.errstate(invalid="ignore"):  # Ignore warnings of division by 0
        CPM = 1e6 * X / sumX

        # Samples with zeros for all genes get nan but should be 0.0
        np.nan_to_num(CPM, copy=False)

    return CPM


def quantile(X):
    """Quantile normalize gene expression to average distribution.

    The procedure is implemented as described on Wikipedia_ and runs on columns and rows:

    * Rearrange the column values so each column is in order from lowest to highest value
    * Find the mean for each row to determine the average distribution of expression values
    * For each column in original data determine a rank from lowest to highest
    * Take the ranking order and substitute in new values from the average distribution

    .. _Wikipedia: https://en.wikipedia.org/wiki/Quantile_normalization
    """
    # Cast array_like objects to Numpy
    X_ = np.asarray(X, dtype=np.float64)
    assert np.min(X_) >= 0.0  # Gene expression must be non-negative

    average_expression_distribution = np.mean(np.sort(X_, axis=0), axis=1)

    rank_avg = rankdata(X_, method="average", axis=0) - 1
    rank_floor = rank_avg.astype(int)
    rank_ceil = np.ceil(rank_avg).astype(int)

    X_floor = average_expression_distribution.take(rank_floor)
    X_ceil = average_expression_distribution.take(rank_ceil)
    X_ = (X_floor + X_ceil) / 2

    if isinstance(X, pd.DataFrame):
        return pd.DataFrame(X_, index=X.index, columns=X.columns)
    else:
        return X_
