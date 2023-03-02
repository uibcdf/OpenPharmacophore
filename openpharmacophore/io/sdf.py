from rdkit import Chem


def load_sdf(file_path):
    """ Load an sdf file with molecules that may contain multiple conformers.

        Parameters
        ----------
        file_path : str

        Returns
        -------
        list : [rdkit.Chem.Mol]
    """
    supp = Chem.SDMolSupplier(file_path, removeHs=False)
    molecules = {}

    for mol in supp:
        name = mol.GetProp("_Name")
        try:
            molecules[name].AddConformer(mol.GetConformer(), assignId=True)
        except KeyError:
            molecules[name] = mol

    return list(molecules.values())
