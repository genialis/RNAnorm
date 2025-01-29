"""Utilities for normalization."""

import numpy as np
from sklearn.preprocessing import FunctionTransformer
from sklearn.utils.validation import validate_data

from ..typing import Numeric2D


class LibrarySize(FunctionTransformer):
    """Library size.

    :param allow_nan: If true, allow X to contain NaN values.

    .. rubric:: Examples

    >>> from rnanorm.datasets import load_toy_data
    >>> from rnanorm import LibrarySize
    >>> X = load_toy_data().exp
    >>> X
              Gene_1  Gene_2  Gene_3  Gene_4  Gene_5
    Sample_1     200     300     500    2000    7000
    Sample_2     400     600    1000    4000   14000
    Sample_3     200     300     500    2000   17000
    Sample_4     200     300     500    2000    2000
    >>> LibrarySize().set_output(transform="pandas").fit_transform(X)
              Library size
    Sample_1       10000.0
    Sample_2       20000.0
    Sample_3       20000.0
    Sample_4        5000.0

    """

    def __init__(self, allow_nan: bool = False) -> None:
        """Initialize class."""
        self.allow_nan = allow_nan
        super().__init__(
            func=np.nansum,
            kw_args=dict(axis=1),
            validate=True,
            feature_names_out=lambda self, input_features: ["Library size"],
        )

    def _check_input(self, X: Numeric2D, *, reset: bool) -> Numeric2D:
        """Check input.

        This function is modified only so that NaN values can be
        tolerated even if self.validate=True.
        """
        kwargs = dict()
        if self.allow_nan:
            kwargs["ensure_all_finite"] = "allow-nan"
        return validate_data(self, X, reset=reset, **kwargs)
