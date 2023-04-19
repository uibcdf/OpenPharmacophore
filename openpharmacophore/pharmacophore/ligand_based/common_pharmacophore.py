from collections import Counter, defaultdict, namedtuple
from dataclasses import dataclass
import heapq
import itertools
from typing import List, Optional, Tuple
import numpy as np
import pyunitwizard as puw

from openpharmacophore import PharmacophoricPoint, Pharmacophore
from openpharmacophore import constants
from openpharmacophore.pharmacophore.align import align_pharmacophores
from openpharmacophore.utils.maths import points_distance, nearest_bins
from openpharmacophore.molecular_systems.chem_feats import ChemFeatContainer


class FeatQueueWrapper:
    """ Wrapper object for a KSublist si it can be stored in a FeatListQueue
    """
    def __init__(self, sublist):
        self.sublist = sublist  # type: KSubList

    def __lt__(self, other):
        return self.sublist.score < other.sublist.score

    def __eq__(self, other):
        return self.sublist.score == other.score


class FeatListQueue:
    """ Class that functions as a min heap for k sublists.
    """

    def __init__(self, size):
        self.size = size
        self._heap = []  # type: list[FeatQueueWrapper]
        self._unique = set()  # type: set[int]

    def reverse_items(self):
        """ Returns all sublist in the queue in descending order of
            score.
        """
        return reversed([it.sublist for it in self._heap])

    def push(self, sublist):
        """ Push a sublist to the queue.

            If the capacity is exceeded the sublist with the smallest score is
            removed and the new sublist inserted. If the sublist to insert has a
            lower score than the smallest element it will not be inserted.

            Only unique sublists are inserted.

            Parameters
            ----------
            sublist : KSubList
        """
        if sublist.id in self._unique:
            return

        if self.size is not None and len(self._heap) >= self.size:
            if sublist.score < self.peek_left().score:
                return
            self.pop()

        heapq.heappush(self._heap, FeatQueueWrapper(sublist))
        self._unique.add(sublist.id)

    def pop(self):
        """ Get the sublists with the smallest score.

            Returns
            -------
            KSubList
        """
        sublist = heapq.heappop(self._heap).sublist
        self._unique.remove(sublist.id)
        return sublist

    def peek_left(self):
        """ Returns the smallest element.

            Returns
            -------
            KSubList
        """
        return self._heap[0].sublist

    def peek_right(self):
        """ Returns the greatest element.

            Returns
            -------
            KSubList
        """
        return self._heap[-1].sublist

    def empty(self):
        """ Check if the queue is empty

            Returns
            -------
            bool
        """
        return len(self._heap) == 0


class ScoringFunction:
    """ A customizable function to score the common pharmacophores.
    """

    def __init__(
            self,
            point_weight=1.0,
            vector_weight=1.0,
            rmsd_cutoff=1.2,  # in angstroms
            cos_cutoff=0.5
    ):
        self.point_weight = point_weight
        self.vector_weight = vector_weight
        self.rmsd_cutoff = rmsd_cutoff
        self.cos_cutoff = cos_cutoff

    def point_score(self, rmsd):
        """ Compute the point score function that is defined by
            1 - RMSD / RMSD_cutoff.

            Where RMSD is the value of the root mean square deviation of the
            alignment between two pharmacophores.

            Parameters
            ----------
            rmsd : float
                RMSD of alignment

            Returns
            -------
            float
        """
        return 1 - rmsd / self.rmsd_cutoff

    def vector_score(self):
        #  TODO: Implement me!
        return 0

    def __call__(self, rmsd):
        return self.point_weight * self.point_score(rmsd) \
            + self.vector_weight * self.vector_score()


class FeatureList:
    """ Class to store feature lists.

       A feature lists represents the pharmacophoric points of a ligand
       including its chemical features and the interpoint distances of its points.

       Attributes
       ----------
       variant : str
           A string with the chemical features of the ligand. For example,
           a conformer with two hydrogen bond acceptors (A) and an aromatic
           ring (R) would have variant 'AAR'.

       distances : np.ndarray
           An array of shape (n_conformers, n_pairs) where n_pairs is the number of interpoint
           distances given by n_points x (n_points - 1) / 2

        coords : QuantityLike
            A quantity with the coordinates of the chemical features.
            Shape (n_conformers, n_feats, 3)
   """

    def __init__(self, variant, distances, coords):
        self.variant = variant
        self.distances = distances
        self.coords = coords

    @classmethod
    def from_chem_feats(cls, chem_feats):
        """ Create a feature list form a molecule chemical features.

            Parameters
            ----------
            chem_feats : list[ChemFeatContainer]
               List of chemical features.

            Returns
            -------
            FeatureList
        """
        variant = ""
        variant += "A" * len(chem_feats[0].acceptor)
        variant += "D" * len(chem_feats[0].donor)
        variant += "H" * len(chem_feats[0].hydrophobic)
        variant += "N" * len(chem_feats[0].negative)
        variant += "P" * len(chem_feats[0].positive)
        variant += "R" * len(chem_feats[0].aromatic)

        n_conformers = len(chem_feats)
        n_feats = len(chem_feats[0])

        coords = puw.quantity(np.zeros((n_conformers, n_feats, 3)), "angstroms")
        for ii in range(len(chem_feats)):
            for jj, feats in enumerate(chem_feats[ii]):
                coords[ii][jj] = feats.coords

        distances = FeatureList.distance_vector(puw.get_value(coords, "angstroms"))

        return cls(variant, distances, coords)

    def has_variant(self, variant):
        """ Returns true if the feature list contains the given variant. """
        counter_self = Counter(self.variant)
        counter_other = Counter(variant)
        for ftype, cnt in counter_other.items():
            cnt_self = counter_self.get(ftype, 0)
            if cnt > cnt_self:
                return False
        return True

    @staticmethod
    def distance_vector(coords):
        """ Compute the vector of inter-site distances.

            Parameters
            ----------
            coords : np.ndarray
               An array with the chemical features positions. Shape
               (n_conformers, n_feats, 3)

            Returns
            -------
            np.array
                An array of shape (n_feats * (n_feats - 1) / 2)
        """
        n_confs = coords.shape[0]
        n_sites = coords.shape[1]
        vec_len = int((n_sites * (n_sites - 1)) / 2)
        distances = np.zeros((n_confs, vec_len))

        for conf in range(n_confs):
            ii = 0
            for site_i in range(n_sites):
                for site_j in range(site_i + 1, n_sites):
                    distances[conf][ii] = points_distance(
                        coords[conf][site_i], coords[conf][site_j]
                    )
                    ii += 1

        return distances

    def subvariant_distances(self, variant, conf):
        """ Get the distance vector of a subvariant.

            Parameters
            ----------
            variant : list[int]
            conf : int

            Returns
            -------
            distances : np.ndarray
        """
        distances = np.zeros(int(len(variant) * (len(variant) - 1) / 2))
        cnt = 0

        n_pairs = self.distances.shape[-1]
        n_feats = len(self.variant)

        for ii in range(len(variant)):
            var_left = variant[ii]
            start = int(n_pairs - (n_feats - var_left) * (n_feats - var_left - 1) / 2)
            for jj in range(ii + 1, len(variant)):
                shift = variant[jj] - var_left - 1
                distances[cnt] = self.distances[conf][start + shift]
                cnt += 1

        return distances

    def k_sublists(self, k_variant):
        """ Get a sublists of the given variant.

            Parameters
            ----------
            k_variant : KVariant
            mol: int

            Returns
            -------
            sublists: list[KSubList]
        """
        n_confs = self.distances.shape[0]
        sublists = [None] * n_confs
        for ii in range(n_confs):
            variant_ind = list(k_variant.var_ind)
            sublists[ii] = KSubList(
                self.subvariant_distances(variant_ind, ii), (k_variant.mol, ii), variant_ind
            )
        return sublists

    def __len__(self):
        return self.coords.shape[0]


@dataclass(frozen=False)
class KSubList:
    """ Class to store sub sets of a feature lists with k points.

        Does not store chem feats coordinates.
    """
    ID = 1

    distances: np.ndarray
    mol_id: Tuple[int, int]
    feat_ind: List[int]
    score: Optional[float] = None

    def __post_init__(self):
        self.id = KSubList.ID
        KSubList.ID += 1

    def __eq__(self, other):
        if other.mol_id != self.mol_id:
            return False
        if other.feat_ind != self.feat_ind:
            return False
        return np.all(self.distances == other.distances)


KVariant = namedtuple("KVariant", ["var_ind", "mol"])


class SurvivingBox:

    def __init__(self):
        self.sublists = []  # type: list[KSubList]
        self.ligands = set()  # type: set[int]
        self.id = []  # type: list[int]

    def __eq__(self, other):
        return self.id == other.id


class CommonPharmacophoreFinder:
    """ Class to search for common pharmacophores in a set of ligands.

        Parameters
        ----------
        scoring_fn_params : dict[str, float], optional
            The parameters of the scoring function.
    """

    def __init__(self, scoring_fn_params=None, **kwargs):
        self.min_dist = kwargs.get("min_dist", puw.quantity(2.0, "angstroms"))
        self.max_dist = kwargs.get("max_dist", puw.quantity(15.0, "angstroms"))
        self.bin_size = kwargs.get("bin_size", puw.quantity(1.0, "angstroms"))

        self.min_dist = puw.get_value(self.min_dist)
        self.max_dist = puw.get_value(self.max_dist)
        self.bin_size = puw.get_value(self.bin_size)

        if scoring_fn_params is None:
            self.scoring_fn = ScoringFunction()
        else:
            self.scoring_fn = ScoringFunction(**scoring_fn_params)

    def find_common_pharmacophores(self, chemical_features, n_points,
                                   min_actives=None, max_pharmacophores=None):
        """ Find common pharmacophores.

            Parameters
            ----------
            chemical_features : list[list[ChemFeatContainer]]
                A nested list where each entry represents the chemical features
                of a ligand and its conformers (as sublist)

            n_points : int
               Extracted pharmacophores will have this number of pharmacophoric
               points.

            min_actives : int, optional
               Number of ligands that must match a common pharmacophore.

            max_pharmacophores : int, optional
               Maximum number of pharmacophores to return. If set to null
               all found pharmacophores will be returned.


            Returns
            -------
            list[Pharmacophore]
                List with the common pharmacophores.

        """
        KSubList.ID = 1

        n_ligands = len(chemical_features)
        if min_actives is None:
            min_actives = n_ligands

        feature_lists = self._get_feat_lists(chemical_features)
        common_variants = self._common_k_point_variants(feature_lists, n_points, min_actives)
        # TODO: eliminate variants with overlaping features
        sub_lists = self._variant_sublists(feature_lists, common_variants)

        scores = {}
        queue = FeatListQueue(size=max_pharmacophores)
        for variant in sub_lists.keys():
            surviving_boxes = self._recursive_partitioning(sub_lists[variant], min_actives)
            for box in surviving_boxes:
                if len(box) > 0:
                    top_representative = self._box_top_representative(box, scores, feature_lists)
                    if top_representative is not None:
                        queue.push(top_representative)

        return self._get_pharmacophores(queue, feature_lists)

    def _recursive_partitioning_util(self, sublists, min_actives, dim, n_dims, boxes):
        """ Helper function for _recursive_partitioning

            Parameters
            ----------
            sublists : list[KSubList]
                Sublists of the same variant

            min_actives : int
                Minimum number of actives that a surviving box must contain

            dim : int
                The dimension of the distances array which will be used for
                partitioning.

            n_dims : int
                Total number of dimensions in the distances array.

            boxes : list[SurvivingBox]
                Stores the surviving sublists
        """
        bins = defaultdict(SurvivingBox)

        for sublist_ in sublists:
            low_bin, upp_bin = nearest_bins(sublist_.distances[dim], self.bin_size)
            if low_bin >= self.max_dist or upp_bin >= self.max_dist:
                continue

            bins[low_bin].sublists.append(sublist_)
            bins[upp_bin].sublists.append(sublist_)

            bins[low_bin].ligands.add(sublist_.mol_id[0])
            bins[upp_bin].ligands.add(sublist_.mol_id[0])

        for box in bins.values():
            if len(box.ligands) >= min_actives:
                if dim < n_dims - 1:
                    self._recursive_partitioning_util(box.sublists, min_actives, dim + 1, n_dims, boxes)
                else:
                    if box not in boxes:
                        boxes.append(box)

    def _recursive_partitioning(self, sublists, min_actives):
        """ Partition feature lists by their inter-site distances and group
            them into "surviving boxes", which are lists of feature lists that
            have similar inter-site distances.

           Parameters
           ----------
           sublists : list[KSubList]
                Sublists of the same variant

           min_actives : int
                Minimum number of actives that a surviving box must contain

           Returns
           -------
           boxes : list[list[KSubList]]
        """
        boxes = []  # type: list[SurvivingBox]
        if sublists:
            n_dims = sublists[0].distances.shape[0]
            self._recursive_partitioning_util(sublists, min_actives, 0, n_dims, boxes)
        return [b.sublists for b in boxes]

    def _box_top_representative(self, box, scores, feature_lists):
        """ Obtain the top ranked representative feature sublist from a surviving box.

            Parameters
            ----------
            box : list[KSublist]
                A surviving box

            scores : dict[tuple, tuple]
                Values of point scores between feature lists.

            Returns
            -------
            KSubList
                The best ranked feature sublist in the box with its score. Can be null if the
                alignments of the feature list are above the rmsd threshold.
        """
        best_score = float("-inf")
        best_ref = None

        for ref in box:
            # Store the scores of the best alignments of the reference against feature lists
            # of other molecules
            best_partner = [float("-inf")] * len(feature_lists)

            exclude = False
            for partner in box:
                # Score reference only with respect to feature lists of other molecules
                ref_mol = ref.mol_id[0]
                partner_mol = partner.mol_id[0]

                if ref_mol == partner_mol:
                    continue

                idx_1 = ref.id
                idx_2 = partner.id
                if idx_2 < idx_1:
                    idx_1, idx_2 = idx_2, idx_1

                try:
                    align_score = scores[(idx_1, idx_2)][0]
                except KeyError:

                    ref_flist, partner_flist = feature_lists[ref_mol], feature_lists[partner_mol]
                    conf_ref, conf_partner = ref.mol_id[1], partner.mol_id[1]

                    align_rmsd, transform = align_pharmacophores(
                        puw.get_value(ref_flist.coords[conf_ref][ref.feat_ind], "angstroms"),
                        puw.get_value(partner_flist.coords[conf_partner][partner.feat_ind], "angstroms"))

                    align_score = self.scoring_fn(align_rmsd)

                    scores[(idx_1, idx_2)] = (align_score, transform)

                if align_score < 0:
                    # Exclude this reference as a potential hypothesis
                    exclude = True
                    break

                best_partner[partner_mol] = max(align_score, best_partner[partner_mol])

            if exclude:
                continue

            # Get average score
            total_score = sum(s for s in best_partner if s > 0) / len([s for s in best_partner if s > 0])
            if total_score > best_score:
                best_score = total_score
                best_ref = ref

        if best_ref is None:
            return
        best_ref.score = best_score
        return best_ref

    @staticmethod
    def _get_pharmacophores(queue, feature_lists):
        """ Extract sublists from the queue and transform them to
            pharmacophores.

            Parameters
            ----------
            queue : FeatListQueue
            feature_lists : list[FeatureList]

            Returns
            -------
            list[Pharmacophore]
        """
        pharmas = []
        radius = puw.quantity(1.0, "angstroms")
        for sublist in queue.reverse_items():
            ref_mol = sublist.mol_id[0]
            ref_conf = sublist.mol_id[1]

            flist = feature_lists[ref_mol]
            coords = flist.coords[ref_conf]

            points = []
            for ii in sublist.feat_ind:
                feat_name = constants.CHAR_TO_FEAT[flist.variant[ii]]
                points.append(PharmacophoricPoint(feat_name, coords[ii], radius))

            pharmas.append(Pharmacophore(points, score=sublist.score, ref_mol=ref_mol, ref_struct=ref_conf))

        return pharmas

    @staticmethod
    def _get_feat_lists(chem_feats):
        """ Get the feat lists of all ligands

            Returns
            -------
            feat_lists : list[FeatureList]
                List where each entry represents the feature lists of
                a ligand
        """
        feat_lists = []
        for cfts in chem_feats:
            feat_lists.append(FeatureList.from_chem_feats(cfts))

        return feat_lists

    @staticmethod
    def _common_k_point_variants(feat_lists, n_points, min_actives):
        """ Find the common k-point variants, that is, the variants
            consisting of k features that are common to at least the
            specified number of actives.

            Parameters
            ----------
            feat_lists : list[FeatureList]
                A feature list for each of the ligands

            n_points : int

            min_actives: int

            Returns
            -------
            k_variants : dict[str, list[KVariant]]

        """
        variant_count = {}
        k_variants = defaultdict(list)

        for ii, flist in enumerate(feat_lists):
            # Keep track of the variants in this ligand
            lig_variants = set()
            for k_var_ind in itertools.combinations(
                    range(len(flist.variant)), n_points):
                k_var = ""
                for jj in k_var_ind:
                    k_var += flist.variant[jj]

                if k_var not in lig_variants:
                    lig_variants.add(k_var)
                    variant_count[k_var] = variant_count.get(k_var, 0) + 1

                k_variants[k_var].append(KVariant(k_var_ind, ii))

        for k_var, count in variant_count.items():
            if count < min_actives:
                k_variants.pop(k_var)

        return k_variants

    @staticmethod
    def _variant_sublists(feat_lists, common_variants):
        """ Group the feature lists by variant type.

            Returns a dictionary of feature lists grouped by variant

            Parameters
            ----------
            feat_lists : list[FeatureList]
                Feature lists of all ligands

            common_variants : dict[str, list[KVariant]]
                The common variants

            Returns
            -------
            variants : dict[str, list[KSubList]]
        """
        var_sublists = defaultdict(list)
        for var, var_list in common_variants.items():
            for k_var in var_list:
                flist = feat_lists[k_var.mol]
                var_sublists[var].extend(flist.k_sublists(k_var))

        return var_sublists

    def __call__(self, chemical_features, n_points, min_actives=None, max_pharmacophores=None):
        """ Find common pharmacophores.

            Shortcut for calling CommonPharmacophoreFinder.find_common_pharmacophores

            Parameters
            ----------
            chemical_features : list[list[ChemFeatContainer]]
        """
        return self.find_common_pharmacophores(chemical_features, n_points, min_actives, max_pharmacophores)
