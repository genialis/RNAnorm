"""Commonly used type hints."""

import typing

import numpy as np
import numpy.typing as npt
from pandera.typing import DataFrame, Series

#: Type for a 1D matrix of numeric values.
Numeric1D = typing.Union[Series[float], npt.NDArray[np.float64]]
#: Type for a 2D matrix of numeric values.
Numeric2D = typing.Union[DataFrame[float], npt.NDArray[np.float64]]

if hasattr(typing, "Self"):
    Self = typing.Self
else:
    Self = object
