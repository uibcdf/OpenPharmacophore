

class DifferentNumAtomsError(ValueError):
    """ Exception raised in fix ligand when the smiles and the
        number of atoms in the ligand don't match.
    """

    def __init__(self, atoms_1, atoms_2):
        self.message = f"Ligand contains {atoms_1} atoms, while " \
                       f"smiles contains {atoms_2} atoms"
        super().__init__(self.message)


class LigandCreationError(ValueError):
    """ Exception raised when failing to create a ligand from a string
    """

    def __init__(self, string, form):
        self.message = f"Failed to create ligand from form {form} with the "\
                        f"string {string}"
        super().__init__(self.message)
