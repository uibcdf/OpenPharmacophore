class NoConformersError(Exception):
    """
        Exception raised when an rdkit molecule has no conformers and thus, no 3D coordinates
    """
    def __init__(self, n_conformers, message="Pharmacophoric points cannot be found without 3D coordinates."):
        self.n_conformers = n_conformers
        self.message = message
        super().__init__(self.message)
    
    def __str__(self):
        return f"Molecule has {self.n_conformers} conformers. -> {self.message}"

class PointTypeError(Exception):
    """
        Exception raised when an invalid pharmacophore point type is passed to a function
    """

    def __init__(self, point_type, message="Invalid point type"):
        self.point_type = point_type
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"{self.message}. \"{self.point_type}\" is not a valid point type"