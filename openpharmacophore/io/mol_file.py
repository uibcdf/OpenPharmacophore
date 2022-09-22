from rdkit import Chem
from openpharmacophore.io import load_mol2_ligands


def mol_file_to_list(file_name):
    """ Load molecules from a file.

        Parameters
        ----------
        file_name: str
            Name of the file.

        Returns
        -------
        list[rdkit.Mol]
            List of ligands.
    """
    ligands = []
    if file_name.endswith(".smi"):
        ligands = [mol for mol in Chem.SmilesMolSupplier(file_name, titleLine=False)]
    elif file_name.endswith(".mol2"):
        ligands = load_mol2_ligands(file_name)
    elif file_name.endswith("sdf"):
        ligands = [mol for mol in Chem.SDMolSupplier(file_name)]
    return ligands
