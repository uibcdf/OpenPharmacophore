from openpharmacophore.io.mol2 import load_mol2_file
from openpharmacophore._private_tools.exceptions import OpenPharmacophoreIOError
from rdkit import Chem
from typing import List

def load_molecules_file(file_name: str, *kwargs) -> List[Chem.Mol]:
    """ Load a file of molecules into a list of rdkit molecules.

        Parameters
        ----------
        file_name : str
            Name of the file that will be loaded.
        
        Returns
        -------
        list of rdkit.Chem.Mol
    """

    if file_name.endswith(".smi"):
        if kwargs:
            if "delimiter" in kwargs and "titleLine" in kwargs:
                ligands = Chem.SmilesMolSupplier(file_name, delimiter=kwargs["delimiter"], titleLine=kwargs["titleLine"])
            elif "delimiter" in kwargs:
                ligands = Chem.SmilesMolSupplier(file_name, delimiter=kwargs["delimiter"], titleLine=True)
            elif "titleLine" in kwargs:
                ligands = Chem.SmilesMolSupplier(file_name, delimiter=' ', titleLine=kwargs["titleLine"])
        else:    
            ligands = Chem.SmilesMolSupplier(file_name, delimiter=' ', titleLine=True)
    elif file_name.endswith(".mol2"):
        ligands = load_mol2_file(file_name)
    elif file_name.endswith(".sdf"):
        ligands = Chem.SDMolSupplier(file_name)
    else:
        raise NotImplementedError
    
    # For some strange reason list(ligands) doesn't work if len(ligands) is not called before
    len(ligands)
    ligands = [lig for lig in ligands if lig is not None]
    if len(ligands) == 0:
        raise OpenPharmacophoreIOError("Molecules couldn't be loaded")
    
    return ligands