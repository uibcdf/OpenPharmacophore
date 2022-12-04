from openpharmacophore.pharmacophore.chem_feats import feature_centroids
from openpharmacophore.utils.maths import points_distance
from openpharmacophore import PharmacophoricPoint
import numpy as np
import pyunitwizard as puw
import math
import itertools
from collections import defaultdict, Counter, namedtuple


BIN_SIZE = 1.0  # In angstroms
MAX_DIST = 15.0  # Maximum interpoint distance in angstroms
MIN_DIST = 2.0  # In angstroms
BINS = np.arange(0, MAX_DIST + BIN_SIZE, step=BIN_SIZE)

# Scoring function values
W_POINT = 1.0
W_VECTOR = 1.0
RMSD_CUTOFF = 1.2
COS_CUTOFF = 0.5


class Ligand:
    """ Class to store ligands used for RDP algorithm.

        Parameters
        ----------
        mol : rdkit.Chem.Mol
    """

    def __init__(self, mol, chem_feats):
        self.mol = mol
        self.feats = chem_feats
        self.variant = ""
        self.feat_count = {}
        self.distances = None

    @property
    def num_confs(self):
        return self.mol.GetNumConformers()

    def _update_variant(self):
        """ Updates the variant, feat_count and distances attributes.
        """
        for feat_type, indices in self.feats.items():
            self.variant += feat_type * len(indices)

        self.variant = "".join(sorted(self.variant, key=str.lower))
        self.feat_count = Counter(self.variant)

        self.distances = np.ones((
            self.mol.GetNumConformers(), len(self.variant), len(self.variant),
        ), dtype=float) * -1

    def interpoint_distances(self, conf):
        """ Calculate the distances between the pharmacophoric points of
            a conformer.

            The distances attribute is updated with the computed distances.

            Parameters
            ----------
            conf : int
                The index of the conformer
        """
        centroids = np.zeros((len(self.variant), 3))
        ii = 0

        feats = "".join(sorted(self.feat_count.keys(), key=str.lower))
        for feat_type in feats:
            for indices in self.feats[feat_type]:
                centroids[ii] = feature_centroids(self.mol, conf, indices)
                ii += 1

        for ii in range(centroids.shape[0]):
            for jj in range(ii + 1, centroids.shape[0]):
                dist = points_distance(centroids[ii], centroids[jj])
                self.distances[conf, ii, jj] = dist
                self.distances[conf, jj, ii] = dist

    def k_distances(self, k_var, conf):
        """ Returns an array with the interpoint distances of the k-variant
            of the specified conformer.

            Parameters
            ----------
            k_var : tuple[int]
                Indices of the k-variant

            conf : int
                Index of the conformer

            Returns
            -------
            dist : np.ndarray
                Array of rank 1 with the interpoint distances of the k-variant.

        """
        if self.distances[conf, 0, 1] == -1:
            self.interpoint_distances(conf)

        shape = (int(len(k_var) * (len(k_var) - 1) / 2), )
        dist = np.zeros(shape)
        ii = 0  # index in dist array
        for jj in range(len(k_var)):
            for kk in range(jj + 1, len(k_var)):
                dist[ii] = self.distances[conf, k_var[jj], k_var[kk]]
                ii += 1
        return dist


class FeatureList:
    """ Class to store feature lists.

        A feature lists represents the pharmacophoric points of a conformer
        including its chemical features and the interpoint distances of its points.

        Parameters
        ----------
        variant : str
            A string with the chemical features of a conformer. For example,
            a conformer with two hydrogen bond acceptors (A) and an aromatic
            ring (R) would have variant 'AAR'.

        var_ind : tuple[int]
            The indices of the features in the conformer.

        fl_id : tuple[int, int]
            The id of this list. The first number is the molecule index and the second
            the conformer index

        distances : np.ndarray
            An array of shape (n_pairs,) where n_pairs is the number of interpoint
            distances given by n_points x (n_points - 1) / 2

        index : int, optional
            An index to use as an identifier for a list
    """
    def __init__(self, variant, var_ind, fl_id, distances, index=None):
        self.id = fl_id
        self.distances = distances
        self.variant = variant
        self.var_ind = var_ind
        self.index = index


class FLContainer:
    """ A container of feature lists of the same variant.
    """
    def __init__(self, n_mols, variant):
        self.variant = variant
        self.n_pairs = int(len(variant) * (len(variant) - 1) / 2)
        # Each sublist "i" stores the feature list of the ith molecule
        self._flists = [[] for _ in range(n_mols)]
        self._num_flists = 0
        self.mols = set()

    def append(self, feat_list):
        """ Append a new feature list to the container.

            Parameters
            ----------
            feat_list : FeatureList
                A feature list with the same variant as the container
        """
        if self.variant == feat_list.variant:
            mol_ind = feat_list.id[0]
            self._flists[mol_ind].append(feat_list)
            self.mols.add(mol_ind)
            self._num_flists += 1
        else:
            raise ValueError("Incorrect variant")

    def append_multiple(self, f_lists):
        """ Append multiple feat lists to the container.

            Parameters
            ----------
            f_lists : list[FeatureList]
        """
        for fl in f_lists:
            self.append(fl)

    def __getitem__(self, item):
        return self._flists[item]

    def __iter__(self):
        self._mol_ind = 0
        self._fl_ind = 0
        return self

    def __next__(self):
        if self._mol_ind < len(self._flists):
            if self._fl_ind < len(self._flists[self._mol_ind]):
                flist = self._flists[self._mol_ind][self._fl_ind]
                self._fl_ind += 1
            else:
                self._mol_ind += 1
                self._fl_ind = 0
                return self.__next__()
            return flist
        else:
            del self._fl_ind
            del self._mol_ind
            raise StopIteration

    def __len__(self):
        return self._num_flists


def nearest_bins(num, bin_size):
    """ Return the index of the nearest bins of the given number
        in the bins array.

        Parameters
        ----------
        num : float
        bin_size : float

        Returns
        -------
        tuple[int, int]
    """
    if num % 1 <= 0.5:
        low_bin = math.floor(num - bin_size)
    else:
        low_bin = math.ceil(num - bin_size)
    return low_bin, low_bin + 1


def recursive_partitioning(container, dim, n_pairs, boxes, min_actives, n_ligs):
    """ Obtain common pharmacophores by recursive distance partitioning.

        Parameters
        ----------
        container : FLContainer
            A container of feature lists of the same variant.
        dim : int
            The dimension of the distances array which will be used for
            partitioning.
        n_pairs : int
            Number of interpoint distance pairs.
        boxes : list[FLContainer]
            A list where surviving boxes will be stored.
        min_actives : int
            Minimum number of actives that a surviving box must contain
        n_ligs : int
            Total number of ligands.

    """
    # TODO: bin attribute in container is not necessary
    bins = [
        FLContainer(n_ligs, container.variant) for _ in range(BINS.shape[0] - 1)
    ]

    for flist in container:
        # Assign each distance to the two closest bins
        low_bin, upp_bin = nearest_bins(flist.distances[dim], BIN_SIZE)
        if low_bin >= len(bins) or upp_bin >= (len(bins)):
            continue
        bins[low_bin].append(flist)
        bins[upp_bin].append(flist)

    for container in bins:
        if len(container.mols) >= min_actives:
            if dim < n_pairs - 1:
                recursive_partitioning(
                    container, dim + 1, n_pairs, boxes, min_actives, n_ligs
                )
            else:
                boxes.append(container)


def pharmacophore_partitioning(container, min_actives, n_ligs):
    """ Partition pharmacophores by their interpoint distances to find common
        pharmacophores.

        Parameters
        ----------
        container : FLContainer
        min_actives : int

        Returns
        -------
        box : list[FLContainer]
    """
    box = []
    recursive_partitioning(container, 0, container.n_pairs, box, min_actives, n_ligs)
    return box


K_VARIANT = namedtuple("K_VARIANT", ["name", "indices", "mol"])


def common_k_point_variants(ligands,  n_points, min_actives):
    """ Find the common k-point feature lists, that is, the lists with variants
        consisting of k pharmacophoric points that are common to at least the
        specified number of actives.

        Parameters
        ----------
        ligands : list[Ligand]
            A list with the ligands.

        n_points : int
            Number of pharmacophoric points

        min_actives : int
            Number of actives that the variant must be present in to be
            considered common

        Returns
        -------
        list[K_VARIANT]
            A list with all the common k-point variants sorted by variant name.

    """
    count = defaultdict(int)
    all_vars = defaultdict(list)

    for ii in range(len(ligands)):
        # Keep track of the variants in this ligand
        mol_variant = {}
        for k_var_indices in itertools.combinations(
                range(len(ligands[ii].variant)), n_points):
            var = ""
            for jj in k_var_indices:
                var += ligands[ii].variant[jj]

            try:
                mol_variant[var] += 1
            except KeyError:
                count[var] += 1
                mol_variant[var] = 1

            variant = K_VARIANT(var, k_var_indices, ii)
            all_vars[var].append(variant)

    common = []
    for var_name, k_vars in all_vars.items():
        if count[var_name] >= min_actives:
            for k_var in k_vars:
                common.append(k_var)

    return common


def common_k_point_feature_lists(ligands, k_variants):
    """ Returns a feature list container for each of the k-point variants.

        Parameters
        ----------
        ligands : list[Ligand]
            The ligands

        k_variants : list[K_VARIANTS]
            A list sorted by variant name

        Returns
        -------
        containers : list[FLContainer]
            List with a container for each variant
    """
    all_containers = []
    if len(k_variants) == 0:
        return all_containers

    prev = k_variants[0].name
    container = FLContainer(len(ligands), prev)
    distances = None
    fl_index = 0

    for k_var in k_variants:
        lig = ligands[k_var.mol]

        if k_var.name != prev:
            all_containers.append(container)
            container = FLContainer(len(ligands), k_var.name)
            prev = k_var.name

        for ii in range(lig.num_confs):
            fl_id = (k_var.mol, ii)
            distances = lig.k_distances(k_var.indices, ii)
            if np.any(distances < MIN_DIST):
                continue

            feat_list = FeatureList(k_var.name, k_var.indices, fl_id, distances)
            feat_list.index = fl_index
            fl_index += 1
            container.append(feat_list)

    if distances is not None and not np.any(distances < MIN_DIST):
        all_containers.append(container)
    return all_containers


def surviving_box_top_representative(surviving_box, point_scores):
    """ Obtain the top ranked representative feature list from a surviving box

        Parameters
        ----------
        surviving_box : FLContainer
            A surviving box

        point_scores : dict[tuple, float]
            Values of point scores between feature lists to avoid repeated
            computations.

        Returns
        -------
        FeatureList
            The best ranked feature list in the box

    """
    pass


def vector_score(id_1, id_2):
    # TODO: complete me!
    return 0


def feat_list_to_pharma(feat_list, ligand):
    """ Retrieve a pharmacophore from a feature list.

        Parameters
        ----------
        feat_list : FeatureList
        ligand : Ligand

        Returns
        -------
        pharmacophore : list[PharmacophoricPoint]
    """
    pharmacophore = []
    for var_ind in feat_list.var_ind:
        feat_name = ligand.variant[var_ind]
        feat_start_ind = ligand.variant.index(feat_name)
        feat_ind = ligand.feats[feat_name][var_ind - feat_start_ind]
        center = feature_centroids(ligand.mol, feat_list.id[1], feat_ind)
        pharmacophore.append(PharmacophoricPoint(
            PharmacophoricPoint.char_to_feature[feat_name],
            puw.quantity(center, "angstroms"),
            puw.quantity(1.0, "angstroms"),
        ))

    return pharmacophore


def find_common_pharmacophores(mols, chem_feats, n_points, min_actives):
    """ Find common pharmacophores in a set of ligands and assigns a score to each one.

        Parameters
        ----------
        mols : list[rdkit.Chem.Mol]
            List of molecules with conformers.

        chem_feats : list[dict]
            List with the chemical features of each molecule

        n_points : int
            Number of pharmacophoric points the common pharmacophores will have.

        min_actives : int
            Number of ligands that the common pharmacophores are present in.

        Returns
        -------
        pharmacophores : list[PharmacophoricPoint]
        scores : list[int]

    """
    ligands = []
    for ii in range(len(mols)):
        lig = Ligand(mols[ii], chem_feats[ii])
        ligands.append(lig)

    k_variants = common_k_point_variants(ligands, n_points, min_actives)
    boxes = common_k_point_feature_lists(ligands, k_variants)

    top_queue = []
    scores = {}
    for box in boxes:
        surviving_boxes = pharmacophore_partitioning(box, min_actives)
        for suv_box in surviving_boxes:
            top_representative = surviving_box_top_representative(suv_box, scores)
        # Check for uniqueness, score and max_pharmacophores
        add_feat_list_to_queue(top_queue, top_representative, max_pharmacophores)

    return [pharmacophore_from_feat_list(t) for t in top_queue]
