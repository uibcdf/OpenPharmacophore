from rdkit import Chem
from rdkit.Chem import AllChem

def generate_conformers(molecule, n_conformers, random_seed=-1, alignment=False):
    """Generate conformers for a molecule
    
        Parameters
        ----------
        molecule: :obj: rdkit.Chem.rdchem.Mol
            Molecule for which conformers will be generated
        
        n_conformers: int
            number of conformers to generate

        random_seed: float or int 
            random seed to use

        alignment: bool
            If true generated conformers will be aligned (Default: False)
        
        Returns
        -------
        molecule: :obj: rdkit.Chem.rdchem.Mol
            Molecule with conformers
    
    """
    molecule = Chem.AddHs(molecule) # Add hydrogens to generate realistic geometries
    cids = AllChem.EmbedMultipleConfs(molecule, numConfs=n_conformers, randomSeed=random_seed)
    
    if alignment:
        AllChem.AlignMolConformers(molecule)
    return molecule