from rdkit import Chem

def smiles_mol_generator(file_object, delimiter=None, header=True, mol_id=True):
    """ A molecule generator for smi and txt files.

        Parameters
        ----------
        file_object : file
            A file object
        
        delimiter : str, default=None
            The delimiter of the files. For smi files it's normally a space.
            If None is passed the delimiter will be any whitespace character.
        
        header : bool
            Whether the file contains a header
            
        mol_id : bool
            Whether the file contains the id or name of the molecules.
        
        Yields
        ------
        mol : rdkit.Chem.Mol
            A molecule.
    """
    if header:
        file_object.readline()
    
    while True:
        line = file_object.readline()
        if not line:
            break
        splitted_line = line.split(delimiter)
        smiles = splitted_line[0]
        mol = Chem.MolFromSmiles(smiles)
        if mol_id:
            mol.SetProp("_Name", splitted_line[1])
        yield mol
        
def smi_has_header_and_id(file_path, delimiter=None):
    """ Function to check if a smi or txt file has a header and/or also contains the ids or name of
        the molecules.
        
        Parameters
        ----------
        file_path : str
            Path to the file
        
        delimiter : str, default=None
            The delimiter of the files. For smi files it's normally a space.
            If None is passed the delimiter will be any whitespace character.
            
        Returns
        -------
        has_header : bool
            Whether the file contains a header
            
        has_id : bool
            Whether the file contains the id or name of the molecules.
    """
    has_id = False
    
    with open(file_path, "r") as fp:
        line = fp.readline().split(delimiter)
        mol = Chem.MolFromSmiles(line[0])
        if mol is None:
            has_header = True
        else:
            has_header = False
            
        line = fp.readline()
        if len(line.split(delimiter)) > 1:
            has_id = True
        
    return has_header, has_id

def mol2_mol_generator(file_object):
    """ A molecule generator for mol2 files.

        Parameters
        ----------
        file_object : file
            A file object
    
        Yields
        ------
        mol : rdkit.Chem.Mol
            A molecule.
    """
    while True:
      
        line = file_object.readline()
        if not line:
            break
        
        if not line.strip():
            continue
        
        # Skip comments
        if line.startswith("#"):
            continue 
        
        if '@<TRIPOS>MOLECULE' not in line:
            mol2_block = "@<TRIPOS>MOLECULE\n"
            mol2_block += line
        else:
            mol2_block = line
        while True:
             line = file_object.readline()
             if '@<TRIPOS>MOLECULE' in line or len(line) == 0:
                 break
             mol2_block += line
        
        yield Chem.MolFromMol2Block(mol2_block)
        