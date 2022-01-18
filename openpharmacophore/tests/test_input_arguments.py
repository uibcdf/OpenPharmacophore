import openpharmacophore._private_tools.exceptions as exc
from openpharmacophore._private_tools.input_arguments import validate_input_quantity, validate_input_array_like
import pyunitwizard as puw
import pytest

def test_validate_input_quantity():
    
    name = "quantity"
    quantity = puw.quantity([2.5, 3.0], unit="angstroms")
    dim = {"[L]": 1}

    with pytest.raises(exc.IsNotQuantityError, match="is not a quantity"):
        validate_input_quantity(1.0, dim, name)
        validate_input_quantity([1.0, 2.0], dim, name)
    
    with pytest.raises(exc.WrongDimensionalityError):
        validate_input_quantity(quantity, {"[T]": 1}, name)

    with pytest.raises(exc.BadShapeError):
        validate_input_quantity(quantity, dim, name, shape=(3,))
    
    with pytest.raises(exc.QuantityDataTypeError):
        validate_input_quantity(quantity, dim, name, dtype=(int))
        validate_input_quantity(puw.quantity(1, unit="angstroms"), dim, name, dtype=(float))


def test_validate_input_array_like():

    name = "array"
    with pytest.raises(exc.NotArrayLikeError):
        validate_input_array_like({1.0}, (1,), name)
    with pytest.raises(exc.BadShapeError):
        validate_input_array_like([1, 2, 3], (4,), name)
        validate_input_array_like((1, 2), (4,), name)



