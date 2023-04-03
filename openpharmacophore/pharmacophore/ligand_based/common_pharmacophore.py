from openpharmacophore.utils.maths import points_distance, nearest_bins
from openpharmacophore import PharmacophoricPoint, Pharmacophore
from openpharmacophore.molecular_systems.chem_feats import ChemFeat, ChemFeatContainer
import numpy as np
import pyunitwizard as puw
import itertools
from collections import defaultdict, Counter, namedtuple


class PriorityQueue:

    def __init__(self, size):
        pass


class ScoringFunction:
    """ A customizable scoring function to score the common pharmacophores.
    """
    def __init__(
            self,
            point_weight=1.0,
            vector_weight=1.0,
            rmsd_cutoff=puw.quantity(1.2, "angstroms"),
            cos_cutoff=0.5
    ):
        self.point_weight = point_weight
        self.vector_weight = vector_weight
        self.rmsd_cutoff = rmsd_cutoff
        self.cos_cutoff = cos_cutoff


class FeatureList:
    """ Class to store feature lists.

       A feature lists represents the pharmacophoric points of a conformer
       including its chemical features and the interpoint distances of its points.

       Attributes
       ----------
       variant : str
           A string with the chemical features of a conformer. For example,
           a conformer with two hydrogen bond acceptors (A) and an aromatic
           ring (R) would have variant 'AAR'.

       mol_id : tuple[int, int]
           The id of this list. The first number is the molecule index and the second
           the conformer index

       distances : np.ndarray
           An array of shape (n_pairs,) where n_pairs is the number of interpoint
           distances given by n_points x (n_points - 1) / 2

        coords : QuantityLike
            A quantity with the coordinates of the chemical features.
   """
    def __init__(self, variant, mol_id, distances, coords):
        self.variant = variant
        self.mol_id = mol_id
        self.distances = distances
        self.coords = coords

    @classmethod
    def from_chem_feats(cls, chem_feats, mol, conf):
        """ Create a feature list form a molecule chemical features.

            Parameters
            ----------
            chem_feats : ChemFeatContainer
               List of chemical features.

            mol : int
                Index of the molecule the chemical features were extracted from.

            conf : int
                Index of the conformer the chemical features were extracted from.

            Returns
            -------
            FeatureList
        """
        variant = ""
        variant += "A" * len(chem_feats.acceptor)
        variant += "D" * len(chem_feats.donor)
        variant += "H" * len(chem_feats.hydrophobic)
        variant += "N" * len(chem_feats.negative)
        variant += "P" * len(chem_feats.positive)
        variant += "R" * len(chem_feats.aromatic)

        n_feats = len(chem_feats)
        coords = puw.quantity(np.zeros((n_feats, 3)), "angstroms")
        for ii, feat in enumerate(chem_feats):
            coords[ii] = feat.coords

        distances = FeatureList.distance_vector(puw.get_value(coords, "angstroms"))

        mol_id = (mol, conf)
        return cls(variant, mol_id, distances, coords)

    @staticmethod
    def distance_vector(coords):
        """ Compute the vector of inter-site distances.

            Parameters
            ----------
            coords : np.ndarray
               An array with the chemical features positions. Shape
               (n_feats, 3)

            Returns
            -------
            np.array
                An array of shape (n_feats * (n_feats - 1) / 2)
        """
        n_sites = coords.shape[0]
        distances = np.zeros(int((n_sites * (n_sites - 1)) / 2))

        ii = 0
        for site_i in range(n_sites):
            for site_j in range(site_i + 1, n_sites):
                distances[ii] = points_distance(coords[site_i], coords[site_j])
                ii += 1

        return distances


class CommonPharmacophoreFinder:
    """ Class to search for common pharmacophores in a set of ligands.

        Parameters
        ----------
        n_points : int
           Extracted pharmacophores will have this number of pharmacophoric
           points.

       min_actives : int, optional
           Number of ligands that must match a common pharmacophore.

       max_pharmacophores : int, optional
           Maximum number of pharmacophores to return. If set to null
           all found pharmacophores will be returned.

        scoring_fn_params : dict[str, float], optional
            The parameters of the scoring function.
    """

    def __init__(
        self, n_points, min_actives=None, max_pharmacophores=None,
        scoring_fn_params=None, **kwargs
    ):
        self.n_points = n_points
        self.min_actives = min_actives
        self.max_pharmacophores = max_pharmacophores

        self.min_dist = kwargs.get("min_dist", puw.quantity(2.0, "angstroms"))
        self.max_dist = kwargs.get("max_dist", puw.quantity(15.0, "angstroms"))
        self.bin_size = kwargs.get("bin_size", puw.quantity(1.0, "angstroms"))
        self.bins = np.arange(
            0, puw.get_value(self.max_dist + self.bin_size), step=puw.get_value(self.bin_size)
        )

        if scoring_fn_params is None:
            self.scoring_fn = ScoringFunction()
        else:
            self.scoring_fn = ScoringFunction(**scoring_fn_params)

    def find_common_pharmacophores(self, chemical_features):
        """ Find common pharmacophores.

            Parameters
            ----------
            chemical_features : list[list[ChemFeatContainer]]
                A nested list where each entry represents the chemical features
                of a ligand and its conformers (as sublist)

            Returns
            -------
            list[Pharmacophore]
                List with the common pharmacophores.

        """
        n_ligands = len(chemical_features)

        feature_lists = self._get_feat_lists(chemical_features)
        common_variants = self._common_k_point_variants(feature_lists)
        sub_lists = self._feat_sub_lists(feature_lists, common_variants)

        scores = {}
        queue = PriorityQueue(size=self.max_pharmacophores)
        for variant in sub_lists:
            surviving_boxes = self._recursive_partitioning(variant, self.min_actives, n_ligands)
            for box in surviving_boxes:
                if box:
                    top_representative = self._box_top_representative(box, scores, n_ligands)
                    queue.push(top_representative)

        return self._get_pharmacophores(queue, chemical_features)

    def _get_feat_lists(self, chem_feats):
        """ Get the feat lists of all conformers
        """
        feat_lists = []
        for ii in range(len(chem_feats)):
            for jj in range(len(chem_feats[ii])):
                feat_lists.append(FeatureList.from_chem_feats(chem_feats[ii][jj], ii, jj))

        return feat_lists

    def __call__(self, chemical_features):
        """ Find common pharmacophores.

            Shortcut for calling CommonPharmacophoreFinder.find_common_pharmacophores

            Parameters
            ----------
            chemical_features : list[list[ChemFeatContainer]]
        """
        return self.find_common_pharmacophores(chemical_features)
