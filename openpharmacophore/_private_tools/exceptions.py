
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


class NotAPharmacophoreError(TypeError):
    """ Exception raised when an object is expected to be a pharmacophore.
    """

    def __init__(self, obj_type=None):

        self.message = ""
        if obj_type:
            self.message += f"Expected a pharmacophore. Got type {obj_type}"

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


class InvalidFileFormat(ValueError):
    """ Exception raised when a file format is not supported or
        is incorrect.
    """

    def __init__(self, file_format):

        self.message = "Invalid file format"
        if file_format:
            self.message += f" {file_format}.\n"
        else:
            self.message += ".\n"

        super().__init__(self.message)


class PDBFetchError(IOError):
    """ Exception raised when a file format is not supported or
        is incorrect.
    """

    def __init__(self, pdb_id, url):

        self.message = f"Error obtaining pdb with id {pdb_id} from {url}"
        super().__init__(self.message)


class ZincDownloadError(IOError):
    """ Exception raised when a file format is not supported or
        is incorrect.
    """

    def __init__(self, status_code, url):

        self.message = f"Error downloading file from {url}.\nStatus code: {status_code}"
        super().__init__(self.message)


class NoLigandError(ValueError):
    """ Exception raised when a protein-ligand complex contains no ligand.
    """

    def __init__(self):

        self.message = f"This molecular system does not contain any ligand."
        super().__init__(self.message)


class NoLigandIndicesError(ValueError):
    """ Exception raised when trying to extract or remove a ligand
        from PLComplex object without obtaining the indices first.
    """

    def __init__(self):
        self.message = f"_lig_indices is empty. Cannot extract/remove."
        super().__init__(self.message)


class SmilesNotFoundError(ValueError):
    """ Exception raised when the smiles of a ligand is not found.
    """

    def __init__(self, ligand):
        self.message = f"Smiles not found for {ligand}"
        super().__init__(self.message)


class DifferentNumAtomsError(ValueError):
    """ Exception raised in fix ligand when the smiles and the
        number of atoms in the ligand don't match.
    """

    def __init__(self, atoms_1, atoms_2):
        self.message = f"Ligand contains {atoms_1} atoms, while " \
                       f"smiles contains {atoms_2} atoms"
        super().__init__(self.message)


class MolGraphError(ValueError):
    """ Exception raised when a rdkit molecule is failed to be created
        from a pdb file.
    """

    def __init__(self, pdb_file):
        self.message = f"Failed to create molecule for {pdb_file}"
        super().__init__(self.message)


class NoConformersError(ValueError):
    """ Exception raised when a molecule contains no conformers.
    """
    pass
