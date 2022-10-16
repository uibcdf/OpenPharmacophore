# Protein ligand interactions
import numpy as np
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import AllChem


# The following values are in nanometers. Taken from PLIP
BS_DIST = 0.85
HYDROPHOBE_DIST = 0.4

chain_names = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K",
               "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V",
               "W", "X", "Y", "Z"]

# Feature definitions taken from Pharmer
smarts_feat_def = {
    '*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]': 'Hydrophobe',
    'N=[CX3](N)-N': 'PosIonizable',
    '[#16!H0]': 'Donor',
    '[#7!H0&!$(N-[SX4](=O)(=O)[CX4](F)(F)F)]': 'Donor',
    '[#7&!$([nX3])&!$([NX3]-*=[!#6])&!$([NX3]-[a])&!$([NX4])&!$(N=C([C,N])N)]': 'Acceptor',
    '[#8!H0&!$([OH][C,S,P]=O)]': 'Donor',
    '[$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])([CH3X4,'
    'CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,'
    'I]': 'Hydrophobe',
    '[$([+,+2,+3])&!$(*[-,-2,-3])]': 'PosIonizable',
    '[$([-,-2,-3])&!$(*[+,+2,+3])]': 'NegIonizable',
    '[$([CH2X4,CH1X3,CH0X2]~[$([!#1]);!$([CH2X4,CH1X3,CH0X2])])]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]': 'Hydrophobe',
    '[$([CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]~[$([CH2X4,CH1X3,CH0X2]~[$([!#1]);!$([CH2X4,CH1X3,CH0X2])])])]~[CH2X4,'
    'CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]': 'Hydrophobe',
    '[$([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(**[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]': 'Hydrophobe',
    '[$([CX3,SX3,PX3](=O)[O-,OH])](=O)[O-,OH]': 'NegIonizable',
    '[$([CX3](=N)(-N)[!N])](=N)-N': 'PosIonizable',
    '[$([NX3]([CX4])([CX4,#1])[CX4,#1])&!$([NX3]-*=[!#6])]': 'PosIonizable',
    '[$([O])&!$([OX2](C)C=O)&!$(*(~a)~a)]': 'Acceptor',
    '[$([SX4,PX4](=O)(=O)[O-,OH])](=O)(=O)[O-,OH]': 'NegIonizable',
    '[$([S]~[#6])&!$(S~[!#6])]': 'Hydrophobe',
    '[C&r3]1~[C&r3]~[C&r3]1': 'Hydrophobe',
    '[C&r4]1~[C&r4]~[C&r4]~[C&r4]1': 'Hydrophobe',
    '[C&r5]1~[C&r5]~[C&r5]~[C&r5]~[C&r5]1': 'Hydrophobe',
    '[C&r6]1~[C&r6]~[C&r6]~[C&r6]~[C&r6]~[C&r6]1': 'Hydrophobe',
    '[C&r7]1~[C&r7]~[C&r7]~[C&r7]~[C&r7]~[C&r7]~[C&r7]1': 'Hydrophobe',
    '[C&r8]1~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]1': 'Hydrophobe',
    '[CH2X4,CH1X3,CH0X2]~[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]': 'Hydrophobe',
    'a1aaaa1': 'Aromatic',
    'a1aaaaa1': 'Aromatic',
    'c1nn[nH1]n1': 'NegIonizable'}


protein_smarts = {
    'C(=O)[O-,OH,OX1]': 'NegIonizable',
    '[#16!H0]': 'Donor',
    '[#7!H0&!$(N-[SX4](=O)(=O)[CX4](F)(F)F)]': 'Donor',
    '[#7&!$([nX3])&!$([NX3]-*=[!#6])&!$([NX3]-[a])&!$([NX4])&!$(N=C([C,N])N)]': 'Acceptor',
    '[#8!H0&!$([OH][C,S,P]=O)]': 'Donor',
    '[$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,'
    'I])([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]([CH3X4,CH2X3,CH1X2,F,Cl,Br,'
    'I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]': 'Hydrophobe',
    '[$(C(N)(N)=N)]': 'PosIonizable',
    '[$([CH2X4,CH1X3,CH0X2]~[$([!#1]);!$([CH2X4,CH1X3,CH0X2])])]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]': 'Hydrophobe',
    '[$([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(**[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]': 'Hydrophobe',
    '[$([O])&!$([OX2](C)C=O)&!$(*(~a)~a)]': 'Acceptor',
    '[$([S]~[#6])&!$(S~[!#6])]': 'Hydrophobe',
    '[$(n1cc[nH]c1)]': 'PosIonizable',
    '[+,+2,+3,+4]': 'PosIonizable',
    '[-,-2,-3,-4]': 'NegIonizable',
    '[CH2X4,CH1X3,CH0X2]~[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]': 'Hydrophobe',
    'a1aaaa1': 'Aromatic',
    'a1aaaaa1': 'Aromatic'
}


def is_ligand_atom(atom):
    """ Check if an atom belongs to a ligand. """
    if not atom.residue.is_water and not atom.residue.is_protein \
            and atom.residue.n_atoms > 4:
        return True
    return False


def find_ligands_in_traj(traj):
    """ Returns a list of ligand ids in a trajectory.

        Parameters
        ----------
        traj : mdtraj.Trajectory

        Returns
        -------
        ligands : list[str]
            A list of the ligands ids in the topology file
    """
    topology = traj.topology

    ligands = []
    for atom in topology.atoms:
        if is_ligand_atom(atom):
            chain = atom.residue.chain.index
            ligand_id = atom.residue.name + ":" + chain_names[chain]
            if ligand_id not in ligands:
                ligands.append(ligand_id)

    return ligands


def get_ligand_atom_indices(traj, ligand_id):
    """ Returns the indices of the atoms of the ligand with the given id
        in the trajectory

        Parameters
        ----------
        traj : mdtraj.Trajectory

        ligand_id : str

        Returns
        -------
        indices : list[int]

    """
    topology = traj.topology
    indices = []

    ligand, chain = ligand_id.split(":")
    chain_index = chain_names.index(chain)
    for atom in topology.atoms:
        if atom.residue.name == ligand and atom.residue.chain.index == chain_index:
            indices.append(atom.index)

    return indices


def ligand_centroid(traj, ligand_id):
    """ Get the centroid of the ligand with the given id."""
    indices = get_ligand_atom_indices(traj, ligand_id)
    coords = traj.xyz[:, indices, :]
    return np.mean(coords, axis=1)[0]


def maximum_distance(centroid, coordinates):
    """ Get the maximum distance from the centroid to the given
        coordinates
    """
    distance = np.sqrt(np.sum(np.power(coordinates - centroid, 2), axis=1))
    return np.amax(distance)


def ligand_maximum_extent(traj, ligand_id):
    """ Computes the maximum extent of the ligand. This is the maximum distance
        from the centroid to any of its atoms.
    """
    centroid = ligand_centroid(traj, ligand_id)
    indices = get_ligand_atom_indices(traj, ligand_id)
    coords = traj.xyz[:, indices, :][0]
    return maximum_distance(centroid, coords)


def get_binding_site_atoms_indices(lig_center, lig_max_extent, coords):
    """ Get the indices of all the atoms that belong to the binding site.

        Parameters
        ----------
        lig_center : np.ndarray of shape (3,)
        lig_max_extent : float
        coords: np.ndarray of shape(n_atoms, 3)

        Returns
        -------
        bs_indices: np.ndarray
    """
    bs_cutoff = lig_max_extent + BS_DIST
    distance = np.sqrt(np.sum(np.power(coords - lig_center, 2), axis=1))
    return np.where(np.logical_and(distance > lig_max_extent, distance <= bs_cutoff))[0]


def chemical_features(molecule, smarts_patterns, bs_indices=None):
    """ Find chemical features in a molecule using smarts patterns.

        Parameters
        ----------
        molecule : rdkit.Mol

        smarts_patterns: dict[str, str]

        bs_indices: np.ndarray, optional

        Returns
        -------
        features dict[str, list[tuple[int]]
            Dictionary with a list of the atom indices that compose each feature.

    """
    features = defaultdict(list)

    # TODO: can we implement this more efficiently?
    for smarts, feat_name in smarts_patterns.items():
        pattern = Chem.MolFromSmarts(smarts)
        atom_indices = molecule.GetSubstructMatches(pattern)
        if len(atom_indices) > 0:
            for indices in atom_indices:
                skip = False
                if bs_indices is not None:
                    for index in indices:
                        if index not in bs_indices:
                            skip = True
                            break
                if not skip:
                    features[feat_name].append(indices)

    return features


def features_centroid(features, coords):
    """ Compute the centroid of the chemical features.

        Parameters
        ----------
        features : dict[str, list[tuple[int]]]

        coords : np.ndarray

        Returns
        -------
        centroids : dict[str, list[np.array]]

    """
    centroids = defaultdict(list)
    for feat, all_indices in features.items():
        for indices in all_indices:
            feat_coords = coords[0, indices, :]
            feat_center = np.mean(feat_coords, axis=0)
            centroids[feat].append(feat_center)

    return centroids


def map_ligand_features_indices(features,  indices_mapper):
    """ Map the indices of the chemical features of the ligand from
        rdkit to mdtraj's.

        Parameters
        ----------
        indices_mapper : dict[int, int]
        features : dict[str, list[tuple]]
    """
    feats_mapped = defaultdict(list)
    for feat_name, all_indices in features.items():
        for indices in all_indices:
            updated_indices = []
            for index in indices:
                updated_indices.append(indices_mapper[index])
            feats_mapped[feat_name].append(tuple(updated_indices))

    return feats_mapped


def fix_bond_order_from_smiles(molecule, smiles):
    """ Assign the correct bond order to a molecule with the corresponding smiles.

        Parameters
        ----------
        molecule : rdkit.Mol
            A molecule with incorrect bond orders.

        smiles: str
            The smiles that is used as a template to fix bond orders

        Returns
        -------
        rdkit.Mol
            Molecule with correct bond orders
    """
    # TODO: if the molecule and the template dont have the same number of
    #  atoms this will fail. We should create a custom exception for that.
    template = AllChem.MolFromSmiles(smiles)
    return AllChem.AssignBondOrdersFromTemplate(template, molecule)


def hydrophobic_pharmacophoric_points(ligand_hyd_centers, protein_hyd_centers):
    """ Returns a list with the centroid of the hydrophobic pharmacophoric points.

        Parameters
        ----------
        ligand_hyd_centers : list[np.array]
            Ligand hydrophobic features centroids.

        protein_hyd_centers : list[np.array]
            Protein hydrophobic features centroids.

        Returns
        -------
        hyd_points : list[np.array]
            Centers of the hydrophobic pharmacophoric points.
    """
    hyd_points = []
    for lig_center in ligand_hyd_centers:
        for protein_center in protein_hyd_centers:
            # If we find two centers whose distance is <= HYDROPHOBE_DIST we add a
            # pharmacophoric point
            if np.sqrt(np.sum(np.power(protein_center - lig_center, 2))) <= HYDROPHOBE_DIST:
                hyd_points.append(lig_center)
                break

    return hyd_points
