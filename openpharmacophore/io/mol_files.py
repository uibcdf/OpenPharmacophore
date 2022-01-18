from openpharmacophore.io.mol2 import load_mol2_file
from openpharmacophore._private_tools.exceptions import OpenPharmacophoreIOError
from rdkit import Chem

def load_molecules_file(file_name, fextension="", **kwargs):
    """ Load a file of molecules into a list of rdkit molecules.

        Parameters
        ----------
        file_name : str
            Name of the file that will be loaded.
        
        fextension : str
            The extension of the file. Must be passed if file argument is a temporary file.
        
        Returns
        -------
        list of rdkit.Chem.mol
    """
    if not fextension:
        fextension = file_name.split(".")[-1]
    
    if fextension == "smi":
        if kwargs:
            if "delimiter" in kwargs and "titleLine" in kwargs:
                ligands = Chem.SmilesMolSupplier(file_name, delimiter=kwargs["delimiter"], titleLine=kwargs["titleLine"])
            elif "delimiter" in kwargs:
                ligands = Chem.SmilesMolSupplier(file_name, delimiter=kwargs["delimiter"], titleLine=True)
            elif "titleLine" in kwargs:
                ligands = Chem.SmilesMolSupplier(file_name, delimiter=' ', titleLine=kwargs["titleLine"])
        else:    
            ligands = Chem.SmilesMolSupplier(file_name, delimiter=' ', titleLine=True)
    elif fextension == "mol2":
        ligands = load_mol2_file(file_name)
    elif fextension == "sdf":
        ligands = Chem.SDMolSupplier(file_name)
    else:
        raise NotImplementedError
    
    # For some strange reason list(ligands) doesn't work if len(ligands) is not called before
    len(ligands)
    ligands = list(ligands)
    ligands = [lig for lig in ligands if lig is not None]
    if len(ligands) == 0:
        raise OpenPharmacophoreIOError("Molecules couldn't be loaded")
    
    return ligands