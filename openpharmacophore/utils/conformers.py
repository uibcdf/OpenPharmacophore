from rdkit.Chem import AddHs
from rdkit.Chem.AllChem import EmbedMultipleConfs


def generate_conformers(molecule, num_conformers=1, random_seed=-1):
    """ Generate conformers for a molecule

       Parameters
       ----------
       molecule : rdkit.Mol
           Molecule for which conformers will be generated.

       num_conformers: int
           Number of conformers to generate

       random_seed : float or int, optional
           Random seed to use.

        Returns
        -------
        rdkit.Mol
            Molecule with conformers. Hydrogens are added as well.

   """
    molecule = AddHs(molecule)
    EmbedMultipleConfs(molecule, numConfs=num_conformers, randomSeed=random_seed)
    return molecule
