from rdkit import Chem
import numpy as np
from dataclasses import dataclass
from typing import Optional, Tuple

from openpharmacophore.constants import QuantityLike, FEAT_TYPES
from openpharmacophore.utils.maths import ring_normal

SMARTS_LIGAND = {
        "aromatic ring": [
            "a1aaaa1",
            "a1aaaaa1"
        ],
        "hydrophobicity": [
            '[$([S]~[#6])&!$(S~[!#6])]',
            '[C&r3]1~[C&r3]~[C&r3]1',
            '[C&r4]1~[C&r4]~[C&r4]~[C&r4]1',
            '[C&r5]1~[C&r5]~[C&r5]~[C&r5]~[C&r5]1',
            '[C&r6]1~[C&r6]~[C&r6]~[C&r6]~[C&r6]~[C&r6]1',
            '[C&r7]1~[C&r7]~[C&r7]~[C&r7]~[C&r7]~[C&r7]~[C&r7]1',
            '[C&r8]1~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]1',
            '[CH2X4,CH1X3,CH0X2]~[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]',
            '*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]',
            '[$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,'
            'I])([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,'
            'CH2X3,CH1X2,F,Cl,Br,I]',
            '[$([CH2X4,CH1X3,CH0X2]~[$([!#1]);!$([CH2X4,CH1X3,CH0X2])])]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]',
            '[$([CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]~[$([CH2X4,CH1X3,CH0X2]~[$([!#1]);!$([CH2X4,CH1X3,'
            'CH0X2])])])]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2] ',
            '[$([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(**[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]',
        ],
        "negative charge": [
            '[$([-,-2,-3])&!$(*[+,+2,+3])]',
            '[$([CX3,SX3,PX3](=O)[O-,OH])](=O)[O-,OH]',
            '[$([SX4,PX4](=O)(=O)[O-,OH])](=O)(=O)[O-,OH]',
            'c1nn[nH1]n1'
        ],
        "positive charge": [
            'N=[CX3](N)-N',
            '[$([+,+2,+3])&!$(*[-,-2,-3])]',
            '[$([CX3](=N)(-N)[!N])](=N)-N',
            '[$([NX3]([CX4])([CX4,#1])[CX4,#1])&!$([NX3]-*=[!#6])]',
        ],
        "hb acceptor": [
            '[#7&!$([nX3])&!$([NX3]-*=[!#6])&!$([NX3]-[a])&!$([NX4])&!$(N=C([C,N])N)]',
            '[$([O])&!$([OX2](C)C=O)&!$(*(~a)~a)]',
        ],
        "hb donor": [
            '[#16!H0]',
            '[#7!H0&!$(N-[SX4](=O)(=O)[CX4](F)(F)F)]',
            '[#8!H0&!$([OH][C,S,P]=O)]',
        ]
    }


SMARTS_PROTEIN = {
        "aromatic ring": [
            "a1aaaa1",
            "a1aaaaa1"
        ],
        "hydrophobicity": [
            '[$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,'
            'I])([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]([CH3X4,CH2X3,CH1X2,F,Cl,Br,'
            'I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]',
            '[$([S]~[#6])&!$(S~[!#6])]',
            '[CH2X4,CH1X3,CH0X2]~[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]',
        ],
        "negative charge": [
            'C(=O)[O-,OH,OX1]',
            '[-,-2,-3,-4]',
        ],
        "positive charge": [
            '[$(C(N)(N)=N)]',
            '[$(n1cc[nH]c1)]',
            '[+,+2,+3,+4]',
        ],
        "hb acceptor": [
            '[#7&!$([nX3])&!$([NX3]-*=[!#6])&!$([NX3]-[a])&!$([NX4])&!$(N=C([C,N])N)]',
            '[$([O])&!$([OX2](C)C=O)&!$(*(~a)~a)]',
        ],
        "hb donor": [
            '[#16!H0]',
            '[#7!H0&!$(N-[SX4](=O)(=O)[CX4](F)(F)F)]',
            '[#8!H0&!$([OH][C,S,P]=O)]',
        ]

}


@dataclass(frozen=True)
class ChemFeat:
    """ Class to store chemical features.

        Attributes
        ----------
        coords: QuantityLike
            3D position of the chemical feature.

        type: str
            Feature type

        atom_indices: tuple[int], optional
            The indices of the atoms that comprise the chemical feature
             in the molecule.

    """
    coords: QuantityLike
    type: str
    atom_indices: Optional[Tuple[int]] = None


@dataclass(frozen=True)
class AromaticRing(ChemFeat):
    type: str = "aromatic ring"
    normal: Optional[np.ndarray] = None


@dataclass(frozen=True)
class HBDonor(ChemFeat):
    type: str = "hb donor"
    hyd: Optional[np.ndarray] = None


class ChemFeatContainer:
    """ A container that stores chemical features by type.
    """

    def __init__(self, feats=None):
        self.aromatic = []
        self.acceptor = []
        self.donor = []
        self.hydrophobic = []
        self.positive = []
        self.negative = []
        self._n_feats = 0

        if feats is not None:
            self.add_feats(feats)

    def add_feats(self, feats):
        """ Add chemical features to the container.

            Parameters
            ----------
            feats : list[ChemFeat]

        """
        for chem_feat in feats:
            feat_list = self._get_feat_list(chem_feat.type)
            feat_list.append(chem_feat)
            self._n_feats += 1

    def has_feat(self, feat_type):
        """ Check if the container has any chemical features of the
            specified type

            Parameters
            ----------
            feat_type : str
                Type of chemical feature.

            Returns
            -------
            bool

        """
        return len(self._get_feat_list(feat_type)) > 0

    def _get_feat_list(self, feat_type):
        if feat_type == "aromatic ring":
            return self.aromatic
        if feat_type == "hb acceptor":
            return self.acceptor
        if feat_type == "hb donor":
            return self.donor
        if feat_type == "hydrophobicity":
            return self.hydrophobic
        if feat_type == "positive charge":
            return self.positive
        if feat_type == "negative charge":
            return self.negative
        raise ValueError(feat_type)

    def __len__(self):
        return self._n_feats

    def __iter__(self):
        # Iterate all features
        yield from self.aromatic
        yield from self.acceptor
        yield from self.donor
        yield from self.hydrophobic
        yield from self.positive
        yield from self.negative


def feature_indices(feat_def, mol):
    """ Get the indices of the atoms that encompass a chemical
        feature.

        Parameters
        ----------
        feat_def : list[str]
            A smart features definition to find chemical features.

        mol : rdkit.Chem.Mol
            Molecule that will be scanned for the desired features.

        Returns
        -------
        list[tuple[int]]
    """
    feat_indices = []

    for smarts in feat_def:
        pattern = Chem.MolFromSmarts(smarts)
        assert pattern is not None, f"{smarts}"
        all_indices = mol.GetSubstructMatches(pattern)
        for indices in all_indices:
            feat_indices.append(indices)

    return feat_indices


def get_indices(mol, feat_def, types=None):
    """ Get indices of different chemical features.

        Parameters
        ----------
        mol : rdkit.Chem.Mol

        feat_def : dict[str, list[str]]
            A list of smarts patterns for each feature

        types : set[str], optional
            The chemical features that will be searched for.

        Returns
        --------
        indices : dict[str, list[tuple[int]]]
            Indices of each chemical feature type.
    """
    if types is None:
        types = FEAT_TYPES

    indices = {}
    for feat_type, smarts_list in feat_def.items():
        if feat_type in types:
            indices[feat_type] = feature_indices(smarts_list, mol)

    return indices


def create_chem_feats(feat_type, indices, coords, store_indices=False):
    """ Returns a list of chemical features of the same type

        Parameters
        ----------
        feat_type : str
            Type of the feature

        indices : list[tuple[int]]
            List where each entry is a tuple with the indices
            of a chemical feature.

        coords : QuantityLike
            A quantity of shape (n_atoms, 3)

         store_indices : bool, default=False
            Whether to store the indices of the atoms in the chemical features.

        Returns
        -------
        list[ChemFeat]

    """
    feats = []
    for ind_tuple in indices:
        feat_coords = np.mean(coords[ind_tuple, :], axis=0)
        if store_indices:
            chem_feat = ChemFeat(coords=feat_coords, type=feat_type, atom_indices=ind_tuple)
        else:
            chem_feat = ChemFeat(coords=feat_coords, type=feat_type)
        feats.append(chem_feat)
    return feats


def mol_chem_feats(indices, coords, store_indices=False):
    """ Get the chemical features of a molecule.

        Parameters
        ----------
        indices : dict[str, list[tuple[int]]]
            Indices of the atoms that encompass the chemical features.

        coords : QuantityLike
            The coordinates of a conformer of the molecule. Has shape
            (n_atoms, 3)

        store_indices : bool, default=False
            Whether to store the indices of the atoms in the chemical features.

        Returns
        -------
        ChemFeatContainer
            A container with chemical features.

    """
    container = ChemFeatContainer()
    for feat, ind_list in indices.items():
        feats = create_chem_feats(feat, ind_list, coords, store_indices=store_indices)
        container.add_feats(feats)

    return container


def donor_chem_feats(donor_indices, coords, mol):
    """ Compute hydrogen bond donor chemical features.

    Parameters
    ----------
    donor_indices : list[tuple[int]]

    coords : QuantityLike

    mol : rdkit.Chem.Mol

    Returns
    -------
    list[HBDonor]

    """
    feats = []
    for ind_tuple in donor_indices:
        centroid = np.mean(coords[ind_tuple, :], axis=0)
        assert len(ind_tuple) == 1  # Donors should be just one atom

        atom = mol.GetAtomWithIdx(ind_tuple[0])
        for nbr in atom.GetNeighbors():
            hyd = coords[nbr.GetIdx()]
            feats.append(
                HBDonor(coords=centroid, hyd=hyd)
            )
    return feats


def aromatic_chem_feats(aromatic_indices, coords):
    """ Compute aromatic chemical features.

    Parameters
    ----------
    aromatic_indices : list[tuple[int]]

    coords : QuantityLike

    Returns
    -------
    feats : list[AromaticRing]

    """
    feats = []
    for ind_tuple in aromatic_indices:
        centroid = np.mean(coords[ind_tuple, :], axis=0)
        normal = ring_normal(ind_tuple, coords, centroid)
        feats.append(
            AromaticRing(coords=centroid, normal=normal)
        )
    return feats
