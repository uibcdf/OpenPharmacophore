from .exceptions import IncorrectShapeError, NotArrayLikeError, WrongDimensionalityError
import numpy as np
import pyunitwizard as puw


def validate_input_array_like(array, shape, name):

    if isinstance(array, np.ndarray):
        if shape != array.shape:
            raise IncorrectShapeError(shape, name)
    elif isinstance(array, (list, tuple, set)):
        if len(shape) > 1 or shape[0] != len(array):
            raise IncorrectShapeError(shape, name)
    else:
        raise NotArrayLikeError(type(array), name)


def validate_input_quantity(quantity, dimensionality, name, shape=None):

    if shape is not None:
        validate_input_array_like(puw.get_value(quantity), shape, name)

    if puw.get_dimensionality(quantity) != dimensionality:
        raise WrongDimensionalityError(dimensionality, name)
