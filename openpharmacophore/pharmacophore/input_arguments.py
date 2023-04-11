from openpharmacophore.pharmacophore.exceptions import IncorrectShapeError, NotArrayLikeError, \
    NotAQuantityError, WrongDimensionalityError
import numpy as np
import pyunitwizard as puw


def validate_input_array_like(array, shape, name=""):
    """ Check whether an array-like object has the correct type, shape or length and data type.

        Parameters
        ----------
        array : obj
            The array-like object

        shape : tuple
            The shape of the array

        name : str, optional
            The name or label of the array-like object

        Raises
        ------
        NotArrayLikeError
            If the object is not array-like.

        IncorrectShapeError
            If the array like object has incorrect shape.

        """
    if isinstance(array, np.ndarray):
        if shape != array.shape:
            raise IncorrectShapeError(shape, name)
    elif isinstance(array, (list, tuple, set)):
        if len(shape) > 1 or shape[0] != len(array):
            raise IncorrectShapeError(shape, name)
    else:
        raise NotArrayLikeError(type(array), name, shape)


def validate_input_quantity(quantity, dimensionality, name="", shape=None):
    """ Check whether a quantity is of the correct dimensionality and shape.

        Parameters
        ----------
        quantity : Quantity
            The quantity that will be validated.

        dimensionality : dict[str, int]
            A dictionary with the expected dimensionality of the quantity.

        name : str, optional
            The name or label of the quantity that is being validated.

        shape :  tuple, optional
            The shape that the quantity is expected to have. If it's a scalar
            should be set to None.

        Raises
        -------
        WrongDimensionalityError
            If the dimensionality doesn't match with the expected one.

        IncorrectShapeError
            If a non-scalar quantity has an incorrect shape.

        NotArrayLike
            If the expected quantity should be of rank > 1 and a scalar is passed.

        """
    if puw.is_quantity(quantity):
        if shape is not None:
            validate_input_array_like(puw.get_value(quantity), shape, name)

        quantity_dim = {
            dim: val for dim, val in puw.get_dimensionality(quantity).items() if val != 0
        }
        if dimensionality != quantity_dim:
            raise WrongDimensionalityError(dimensionality, name)
    else:
        raise NotAQuantityError(type(quantity), name)
