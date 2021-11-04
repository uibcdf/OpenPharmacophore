from rdkit import Chem
from rdkit.Chem import AllChem

from openpharmacophore._private_tools.exceptions import NoConformersError

def generate_conformers(molecule, n_conformers, random_seed=-1, alignment=False):
    """Generate conformers for a molecule
    
        Parameters
        ----------
        molecule : rdkit.Chem.Mol
            Molecule for which conformers will be generated.
        
        n_conformers: int
            Number of conformers to generate

        random_seed : float or int, optional 
            Random seed to use.

        alignment : bool, optional
            If true generated conformers will be aligned (Default: False).
        
        Returns
        -------
        molecule : rdkit.Chem.Mol
            Molecule with conformers.
    
    """
    if not isinstance(n_conformers, int):
        raise TypeError("n_conformers must be an integer")
    if n_conformers < 0:
        raise ValueError("n_conformers must be greater than 0")
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
        molecule : rdkit.Chem.mol
            The molecule which energy will be calculated.

        forcefield : {"UFF", "MMFF"}, optional.
            The forcefield that will be used to calculate the energy 
            (default="UFF").

        Returns
        -------
        double
            The energy of the molecule.
    """
    if molecule.GetNumConformers() == 0:
        raise NoConformersError("Molecule must have at least one conformer")

    if forcefield == "UFF":
        ff = AllChem.UFFGetMoleculeForceField(molecule, confId=conformer_id)
    elif forcefield == "MMFF":
        props = AllChem.MMFFGetMoleculeProperties(molecule)
        ff = AllChem.MMFFGetMoleculeForceField(molecule, props, confId=conformer_id)

    return ff.CalcEnergy()