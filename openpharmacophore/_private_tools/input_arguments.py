from openpharmacophore._private_tools.exceptions import (
    BadShapeError, QuantityDataTypeError, WrongDimensionalityError, 
    IsNotQuantityError, NotArrayLikeError)
import numpy as np
import pint
import pyunitwizard as puw


def validate_input_quantity(quantity, dimensionality, name, 
                            dtype=(float, int, np.int64, np.float32), 
                            shape=None):
    """ Check whether a quantity is of the correct dimenstionality, shape and data type

        Parameters
        ----------
        quantity : pyunitwizard.Quantity
            The quantity object that will be validated.
        
        dimensionality : dict
            A dictionary with the expected dimensionality of the quantity.
        
        name : str
            The name or label of the quantity that is being validated.
            
        dtype : tuple, default=(float, int, np.int64, np.float32)
            The type or types that the quantity is expected to have
        
        shape :  tuple, optional
            The shape that the quantity is expected to have. If the quantity is a scalar
            this argument shouldn't be passed.
            
        Raises
        -------
        IsNotQuantityError
            If it is not a pyunitwizar quantity.
        
        WrongDimensionalityError
            If the dimensionality doesn't match with the required one.
        
        QuantityDataTypeError
            If the datatype is not of the expected types.
        
        BadShapeError
            If a non scalar quantity has an incorrect shape.
        
    """
    if not isinstance(quantity, pint.Quantity):
        raise IsNotQuantityError(f"{name} is not a quantity")    
    
    quantity_dim = {dim: val for dim, val in puw.get_dimensionality(quantity).items() if val != 0}

    if quantity_dim != dimensionality:
        raise WrongDimensionalityError(f"{name} has dimensionality {quantity_dim}. Expected {dimensionality}.")

    quantity_val = puw.get_value(quantity)
    if isinstance(quantity_val, np.ndarray):
        quantity_val = quantity_val[0]
    if not isinstance(quantity_val, dtype):
        raise QuantityDataTypeError(f"{name} has data type {type(quantity_val)}. Expected {dtype}.")
    
    if shape and isinstance(puw.get_value(quantity), np.ndarray):
        quantity_shape = puw.get_value(quantity).shape
        if quantity_shape != shape:
            raise BadShapeError(f"{name} has shape {quantity_shape}. Expected Shape {shape}")
        
def validate_input_array_like(array, shape, name):
    """ Check whether an array-like object has the correct type, shape or length and data type.
    
        Parameters
        ----------
        array : obj
            The array-like object
        
        shape : tuple
            The shape or lenght of the array
        
        name : str
            The name or label of the array-like object
        
        types : tuple, default=(list, tuple, np.ndarray)
            The valid array-like types.
        
        Raises
        ------
        NotArrayLikeError
            If the object is not array like.
        
        BadShapeError
            If the array like object has incorrect shape/lenght.
    
    """
    types=(list, tuple, np.ndarray)

    if not isinstance(array, types):
        raise NotArrayLikeError(f"{name} is not array like.")
    
    if isinstance(array, np.ndarray):
        array_shape = array.shape
        if array_shape != shape:
            raise BadShapeError(f"{name} has shape {array_shape}. Expected Shape {shape}")
    elif isinstance(array, (list, tuple)):
        array_length = len(array)
        if array_length != shape[0]:
            raise BadShapeError(f"{name} has shape {array_length}. Expected Shape {shape}")
        