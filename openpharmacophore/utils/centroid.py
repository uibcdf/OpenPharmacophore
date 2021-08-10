import numpy as np

def feature_centroid(molecule, atom_indxs, conformer_idx):
    """
        Get the 3D coordinates of the centroid of a feature that encompasses more than 
        one atom. This could be aromatic, hydrophobic, negative and positive features

        Parameters
        ----------
        molecule: :obj: rdkit.Chem.rdchem.Mol
                Molecule that contains the feature which centroid will be computed

        atom_indxs: tuple of int
                Indices of the atoms that belong to the feature

        conformer_idx: int 
                Index of the conformer for which the feature centroid will be computed

        Returns
        -------
        centroid: numpy.ndarray; shape: (3, )
            3D coordinates of the centroid of the feature.

    """
    
    n_atoms = len(atom_indxs)
    coords = np.zeros((n_atoms, 3))
    for j, idx in enumerate(atom_indxs):
            position = molecule.GetConformer(conformer_idx).GetAtomPosition(idx)
            coords[j, 0] = position.x
            coords[j, 1] = position.y
            coords[j, 2] = position.z
    
    centroid = coords.mean(axis=0)
   
    return centroid