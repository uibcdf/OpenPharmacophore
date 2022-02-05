from rdkit import Chem
from typing import List

def load_mol2_file(file_name: str) -> List[Chem.Mol]:
    """ Load molecules from a mol2 file.

        Parameters
        ----------
        file_name : str
            Name of the file containing the ligands
        
        Returns
        ---------
        molecules : list of rdkit.Chem.Mol
       
    """
    molecules = []
    with open(file_name, 'r') as f:
        doc = [line for line in f.readlines()]

    start = [index for (index,p) in enumerate(doc) if '@<TRIPOS>MOLECULE' in p]
    finish = [index-1 for (index,p) in enumerate(doc) if '@<TRIPOS>MOLECULE' in p]
    finish.append(len(doc))
    
    interval = list(zip(start, finish[1:]))
    
    for i in interval:
        block = ",".join(doc[i[0]:i[1]]).replace(',','')
        mol = Chem.MolFromMol2Block(block)
        molecules.append(mol)
    
    return molecules