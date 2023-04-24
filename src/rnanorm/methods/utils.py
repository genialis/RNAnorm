"""Utils for normalization"""
import numpy as np
from sklearn.preprocessing import FunctionTransformer

from ..typing import Numeric2D


class LibrarySize(FunctionTransformer):
    """Get library size of an expression matrix.

    Samples are in rows, genes in columns.

    :param validate: If true, convert X to Numpy array and make
        various checks (2D matrix, non-empty, all values are finite)
        NaN values are not allowed.

    """

    def __init__(self, allow_nan: bool = False) -> None:
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
