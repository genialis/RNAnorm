"""Utilities for normalization."""
import numpy as np
from sklearn.preprocessing import FunctionTransformer

from ..typing import Numeric2D


class LibrarySize(FunctionTransformer):
    """Library size.

    :param allow_nan: If true, allow X to contain NaN values.

    .. rubric:: Examples

    >>> from rnanorm.datasets import load_rnaseq_toy
    >>> from rnanorm import LibrarySize
    >>> X = load_rnaseq_toy(as_frame=True).exp
    >>> X
           G1     G2      G3      G4       G5
    S1  200.0  300.0   500.0  2000.0   7000.0
    S2  400.0  600.0  1000.0  4000.0  14000.0
    S3  200.0  300.0   500.0  2000.0  17000.0
    S4  200.0  300.0   500.0  2000.0   2000.0
    >>> LibrarySize().set_output(transform="pandas").fit_transform(X)
        Library size
    S1       10000.0
    S2       20000.0
    S3       20000.0
    S4        5000.0

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
            kwargs["force_all_finite"] = "allow-nan"
        return self._validate_data(X, reset=reset, **kwargs)
