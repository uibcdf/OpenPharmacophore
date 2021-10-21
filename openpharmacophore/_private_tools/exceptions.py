
class ShapeWithNoColorError(ValueError):

    def __init__(self, documentation_web=None):

        message = ('Either a feture name with a color palette or a color is needed to show'
        'this pharmacophoric shape. Check the online documentation for more information.')

        if documentation_web is not None:
            message = message[:-1] + ': {}'.format(documentation_web)

        super().__init__(message)

class NoConformersError(Exception):
    """
        Exception raised when an rdkit molecule has no conformers and thus, no 3D coordinates.
    """
    def __init__(self, n_conformers, message="Pharmacophoric points cannot be found without 3D coordinates."):
        self.n_conformers = n_conformers
        self.message = message
        super().__init__(self.message)
    
    def __str__(self):
        return f"Molecule has {self.n_conformers} conformers. -> {self.message}"

class InvalidFeatureError(Exception):
    """
        Exception raised when trying to remove an invalid feature from a pharmacophore object.
    """

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

class InvalidFileError(Exception):
    """
        Exception raised when trying to load or save from or to an invalid file type.
    """

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

class MissingParameters(Exception):
    """
        Exception raised when calling a function and not passing all required parameters.
    """

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

class InvalidFileFormat(Exception):
    """
        Exception raised when passing an invalid file format to a function that requires it.
    """

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

class FetchError(Exception):
    """
        Exception raised when fetching a file from a url.
    """

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

