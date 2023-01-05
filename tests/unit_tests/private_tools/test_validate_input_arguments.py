from openpharmacophore.point.input_arguments import validate_input_quantity, validate_input_array_like
from openpharmacophore._private_tools import exceptions as exc
import numpy as np
import pyunitwizard as puw
import pytest


@pytest.fixture()
def scalar_quantity():
    return puw.quantity(1.0, "angstroms")


@pytest.fixture()
def array_quantity():
    return puw.quantity([1.0, 1.0], "angstroms")


def test_validate_input_quantity_wrong_dimensionality_scalar(scalar_quantity):
    validate_input_quantity(scalar_quantity, {"[L]": 1}, "radius")  # should not raise

    with pytest.raises(exc.WrongDimensionalityError,
                       match="radius has incorrect dimensionality"):
        validate_input_quantity(scalar_quantity, {"[L]": 2}, "radius")


def test_validate_input_quantity_wrong_dimensionality_array(array_quantity):
    with pytest.raises(exc.WrongDimensionalityError,
                       match="center has incorrect dimensionality"):
        validate_input_quantity(array_quantity, {"[L]": 2}, "center")


def test_validate_input_quantity_wrong_shape_scalar(scalar_quantity):
    validate_input_quantity(scalar_quantity, {"[L]": 1}, "radius", shape=None)  # Should not raise

    # scalar should raise not array like error
    with pytest.raises(exc.NotArrayLikeError,
                       match="radius is of type <class 'float'>"):
        validate_input_quantity(scalar_quantity, {"[L]": 1}, "radius", shape=(2,))


def test_validate_input_quantity_wrong_shape_array(array_quantity):
    with pytest.raises(exc.IncorrectShapeError,
                       match="center has incorrect shape"):
        validate_input_quantity(array_quantity, {"[L]": 1}, "center", shape=(2, 2))


def test_validate_input_quantity_object_is_not_quantity():
    with pytest.raises(exc.NotAQuantityError, match="Expected a quantity"):
        validate_input_quantity(np.array([1, 2]), {"[L]": 1}, shape=(2,))


def test_validate_input_array_like_with_list():
    with pytest.raises(exc.IncorrectShapeError,
                       match="center has incorrect shape"):
        validate_input_array_like([1, 2], (3,), "center")


def test_validate_input_array_like_with_numpy_array():
    validate_input_array_like(np.array([1, 2]), (2,), "array")  # should not raise

    with pytest.raises(exc.IncorrectShapeError,
                       match="array has incorrect shape"):
        validate_input_array_like(np.array([1, 2]), (3,2), "array")
