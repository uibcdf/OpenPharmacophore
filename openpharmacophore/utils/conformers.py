from rdkit import Chem
from rdkit.Chem import AllChem

from openpharmacophore._private_tools.exceptions import NoConformersError

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

def conformer_energy(molecule, conformer_id=0, forcefield="UFF"):
    """ Get the energy of a conformer in a molecule using the universal forcefield (UFF) forcefield
        or the merck molecular forcefield (MMFF) forcefield.

        Parameters
        ----------
        molecule: rdkit.Chem.mol
            The molecule which energy will be calculated

        forcefield: str, optional, default "UFF".
            The forcefield that will be used

        Returns
        -------
        double
            The energy of the molecule
    """
    if molecule.GetNumConformers() == 0:
        raise NoConformersError("Molecule must have at least one conformer")

    if forcefield == "UFF":
        ff = AllChem.UFFGetMoleculeForceField(molecule, confId=conformer_id)
    elif forcefield == "MMFF":
        props = AllChem.MMFFGetMoleculeProperties(molecule)
        ff = AllChem.MMFFGetMoleculeForceField(molecule, props, confId=conformer_id)

    return ff.CalcEnergy()