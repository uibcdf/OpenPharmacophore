
class WrongDimensionalityError(ValueError):
    """ Exception raised when a quantity does not have
        the expected dimensionality.
    """
    def __init__(self, dimensionality=None, name=""):
        if name:
            self.message = f"{name} has incorrect dimensionality.\n"
        else:
            self.message = f"Incorrect dimensionality.\n"

        if dimensionality is not None:
            self.message += f"Expected {dimensionality}"

        super().__init__(self.message)


class IncorrectShapeError(ValueError):
    """ Exception raised when an array or array like object does not have
        the expected shape.
    """

    def __init__(self, shape=None, name=""):
        if name:
            self.message = f"{name} has incorrect shape.\n"
        else:
            self.message = f"Incorrect shape.\n"

        if shape is not None:
            self.message += f"Expected {shape}"

        super().__init__(self.message)


class NotArrayLikeError(TypeError):
    """ Exception raised when an object is expected to be a numpy
        array, list, tuple or set.
   """

    def __init__(self, obj_type=None, name="", shape=None):

        if shape is not None:
            self.message = f"Expected array-like of shape {shape}.\n"
        else:
            self.message = "Expected array-like.\n"

        if name and obj_type:
            self.message += f"{name} is of type {obj_type}."
        elif obj_type:
            self.message += f"Got type {obj_type}"

        super().__init__(self.message)


class NotAQuantityError(TypeError):
    """ Exception raised when an object is expected to be a numpy
        array, list, tuple or set.
   """

    def __init__(self, obj_type=None, name="", dimensionality=None):

        if dimensionality is not None:
            self.message = f"Expected a quantity of dimensionality {dimensionality}.\n"
        else:
            self.message = "Expected a quantity.\n"

        if name and obj_type:
            self.message += f"{name} is of type {obj_type}."
        elif obj_type:
            self.message += f"Got type {obj_type}"

        super().__init__(self.message)


class InvalidFeatureError(ValueError):
    """ Exception raised when a pharmacophoric feature is not supported.
    """

    def __init__(self, feat_name):

        self.message = "Invalid feature name"
        if feat_name:
            self.message += f" {feat_name}.\n"
        else:
            self.message += ".\n"

        super().__init__(self.message)
