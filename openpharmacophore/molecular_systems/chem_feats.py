from rdkit import Chem
from dataclasses import dataclass
import numpy as np

from openpharmacophore.config import QuantityLike


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
        ]
    }


@dataclass
class ChemFeat:
    type: str
    coords: QuantityLike


class ChemFeatContainer:

    def __init__(self):
        self.aromatic = []
        self.acceptor = []
        self.donor = []
        self.hydrophobic = []
        self.positive = []
        self.negative = []

    def add_feats(self, feats):
        """ Add chemical features to the container.

            Parameters
            ----------
            feats : list[ChemFeat]

        """
        for chem_feat in feats:
            feat_list = self._get_feat_list(chem_feat.type)
            feat_list.append(chem_feat)

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


def feature_centroids(mol, conf, indices):
    """ Get the centroid of a chemical feature.
        Parameters
        ----------
        mol : rdkit.Chem.Mol
            A molecule
        conf : int
            Index of conformer
        indices : tuple[int]
            Atom indices of the chemical feature.
        Returns
        -------
        np.ndarray
            Array of shape (3,)
    """
    coords = np.zeros((len(indices), 3))
    conformer = mol.GetConformer(conf)
    for ii, atom in enumerate(indices):
        pos = conformer.GetAtomPosition(atom)
        coords[ii] = pos.x, pos.y, pos.z

    return np.mean(coords, axis=0)


def create_chem_feats(feat_type, indices, coords):
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

        Returns
        -------
        list[ChemFeat]

    """
    feats = []
    for ind_tuple in indices:
        feat_coords = np.mean(coords[ind_tuple, :], axis=0)
        feats.append(
            ChemFeat(feat_type, feat_coords)
        )
    return feats


def mol_chem_feats(mol, coords, feat_def):
    """ Get the chemical features of a molecule.

        Parameters
        ----------
        mol : rdkit.Chem.Mol
            Molecule which chemical features will be extracted.

        coords : QuantityLike
            The coordinates of a conformer of the molecule. Has shape
            (n_atoms, 3)

        feat_def : {"protein", "ligand"}
            Feature definition that is used to get chemical features.

        Returns
        -------
        ChemFeatContainer
            A container with chemical features.

    """
    if feat_def == "protein":
        smarts_def = SMARTS_PROTEIN
    elif feat_def == "ligand":
        smarts_def = SMARTS_LIGAND
    else:
        raise ValueError(feat_def)

    container = ChemFeatContainer()
    for feat, smarts in smarts_def.items():
        indices = feature_indices(smarts, mol)
        feats = create_chem_feats(feat, indices, coords)
        container.add_feats(feats)

    return container
