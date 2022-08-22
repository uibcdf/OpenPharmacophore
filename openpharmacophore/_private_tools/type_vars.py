# This module will store variables that will serve for type hints
import typing
import numpy as np
from pint import Quantity as pint_quantity
from openmm.unit.quantity import Quantity as openmm_quantity

# We import some of the most common classes of the typing
# library to avoid importing it in each module.
List = typing.List
Optional = typing.Optional
Sequence = typing.Sequence

# Custom types
ArrayLike = typing.TypeVar("ArrayLike",
                           list,
                           np.ndarray)

Quantity = typing.TypeVar("Quantity",
                          pint_quantity,
                          openmm_quantity)
