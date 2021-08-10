import numpy as np

def donor_acceptor_direction_vector(molecule, atom_indx, coords, conformer_idx):
    """
        Compute the direction vector for an H bond donor or H bond acceptor feature 

        Parameters
        ----------
        molecule: :obj: rdkit.Chem.rdchem.Mol
                Molecule that contains the feature which direction vector will be computed

        atom_indx: int
                Index of the H bond acceptor or donor atom

        coords: numpy.ndarray; shape(3,)
                Coordiantes of the H bond acceptor or donor atom

        conformer_idx: int 
                Index of the conformer for which the direction vector will be computed

        Returns
        -------
        direction: numpy.ndarray; shape(3,)
                Coordinates of the direction vector

    """
    direction = np.zeros((3,)) 
    atom = molecule.GetAtomWithIdx(atom_indx)
    for a in atom.GetNeighbors():
        if a.GetSymbol() == "H":
            continue
        position = molecule.GetConformer(conformer_idx).GetAtomPosition(a.GetIdx())
        direction[0] += position.x - coords[0]
        direction[1] += position.y - coords[1]
        direction[2] += position.z - coords[2]
    
    direction = -direction
    return direction
    

def aromatic_direction_vector(molecule, atom_indxs, conformer_idx):
    """
        Compute the direction vector for an aromatic feature 

        Parameters
        ----------
        molecule: :obj: rdkit.Chem.rdchem.Mol
                Molecule that contains the feature which direction vector will be computed

        atom_indx: int
                Index of the H bond acceptor or donor atom


        conformer_idx: int 
                Index of the conformer for which the direction vector will be computed

        Returns
        -------
        direction: numpy.ndarray; shape(3,)
                Coordinates of the direction vector

    """
    coords = np.zeros((3, 3)) # Take just the first three atoms
    for j, idx in enumerate(atom_indxs[0:3]):
            position = molecule.GetConformer(conformer_idx).GetAtomPosition(idx)
            coords[j, 0] = position.x
            coords[j, 1] = position.y
            coords[j, 2] = position.z
    
    # Find the vector normal to the plane defined by the three atoms
    u = coords[1, :] - coords[0, :]
    v = coords[2, :] - coords[0, :]
    direction = np.cross(u, v)

    return direction