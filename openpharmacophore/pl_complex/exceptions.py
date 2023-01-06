

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


class NoLigandError(ValueError):
    """ Exception raised when a protein-ligand complex contains no ligand.
    """

    def __init__(self):

        self.message = f"This molecular system does not contain any ligand."
        super().__init__(self.message)