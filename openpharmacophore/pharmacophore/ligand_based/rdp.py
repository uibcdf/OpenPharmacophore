from openpharmacophore.pharmacophore.chem_feats import feature_centroids
from openpharmacophore.utils.maths import points_distance
from openpharmacophore import PharmacophoricPoint
import numpy as np
import pyunitwizard as puw
import math
import itertools
from collections import defaultdict, Counter, namedtuple


BIN_SIZE = 1.0  # In angstroms
# Maximum interpoint distance in angstroms
MAX_DIST = 15.0
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

        fl_id : tuple[int, int]
            The id of this list. The first number is the molecule index and the second
            the conformer index

        distances : np.ndarray
            An array of shape (n_pairs,) where n_pairs is the number of interpoint
            distances given by n_points x (n_points - 1) / 2
    """
    count = 0  # keep track of how many objects there are

    def __init__(self, variant, var_ind, fl_id, distances, index=None):
        self.n_pairs = len(variant) * (len(variant) - 1) / 2
        if distances.shape[0] != self.n_pairs:
            raise ValueError

        self.id = fl_id
        self.distances = distances
        self.variant = variant
        self.var_ind = var_ind
        self.index = index


class FLContainer:
    """ A container of feature lists of the same variant."""
    def __init__(self, variant="", bin=None):
        self.variant = variant
        self.bin = bin
        self.flists = []
        # Stores the ids of the mols
        self.mols = set()

    def append(self, feat_list):
        """ Append a new feature list to the container.

            Parameters
            ----------
            feat_list : FeatureList
                A feature list with the same variant as the container
        """
        if self.variant == "":
            self.variant = feat_list.variant

        if self.variant == feat_list.variant:
            self.flists.append(feat_list)
            self.mols.add(feat_list.id[0])
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
        return self.flists[item]

    def __iter__(self):
        self.ii = 0
        return self

    def __next__(self):
        if self.ii < len(self.flists):
            flist = self.flists[self.ii]
            self.ii += 1
            return flist
        else:
            raise StopIteration

    def __len__(self):
        return len(self.flists)


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


def recursive_partitioning(container, dim, n_pairs, boxes, min_actives):
    """ Obtain common pharmacophores by recursive distance partitioning.

        Parameters
        ----------
        container : FLContainer
        dim : int
        n_pairs : int
        boxes : list[FLContainer]
        min_actives : int

    """
    # TODO: bin attribute in container is not necessary
    bins = [
        FLContainer(bin=(BINS[ii], BINS[ii + 1])) for ii in range(BINS.shape[0] - 1)
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
                    container, dim + 1, n_pairs, boxes, min_actives
                )
            else:
                boxes.append(container)


def pharmacophore_partitioning(container, min_actives):
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
    recursive_partitioning(container, 0, container[0].n_pairs, box, min_actives)
    return box


def score_common_pharmacophores(box, rmsd_dict):
    """ Calculate the score of all feature lists in a container and sort them by score.

        Parameters
        ----------
        box : FLContainer
            A surviving box

        rmsd_dict : dict[tuple, float]
            Values of rmsd between feature lists to avoid repeated
            computations.

        Returns
        -------
        scores [list[tuple]]
            A list of tuples with the scores between each feature list pairs, ordered in descending
            order of score.

    """
    scores = []
    exclude = []
    for ii in range(len(box)):
        for jj in range(ii + 1, len(box)):
            try:
                rmsd = rmsd_dict[(box[ii].index, box[jj].index)]
            except KeyError:
                rmsd = np.sqrt(((box[ii].distances - box[jj].distances) ** 2).mean())
                rmsd_dict[(box[ii].index, box[jj].index)] = rmsd

            if rmsd >= RMSD_CUTOFF:
                exclude.append(ii)
                break
            point_score = 1 - rmsd / RMSD_CUTOFF
            # vec_score = vector_score(box[ii].id, box[jj].id)
            # score = W_POINT * point_score + W_VECTOR * vec_score
            scores.append((point_score, ii, jj))

    scores.sort(reverse=True)
    if len(exclude) > 0:
        # Filter pharmacophores that exceed rmsd cutoff
        return [s for s in scores if s[1] not in exclude]
    else:
        return scores


def vector_score(id_1, id_2):
    # TODO: complete me!
    return 0


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
    container = FLContainer(prev)
    for k_var in k_variants:
        lig = ligands[k_var.mol]

        if k_var.name != prev:
            all_containers.append(container)
            container = FLContainer(k_var.name)
            prev = k_var.name

        for ii in range(lig.num_confs):
            fl_id = (k_var.mol, ii)
            distances = lig.k_distances(k_var.indices, ii)
            feat_list = FeatureList(k_var.name, k_var.indices, fl_id, distances)
            feat_list.index = FeatureList.count
            FeatureList.count += 1
            container.append(feat_list)

    FeatureList.count = 0
    all_containers.append(container)
    return all_containers


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
    assert len(mols) == len(chem_feats)

    ligands = []
    for ii in range(len(mols)):
        lig = Ligand(mols[ii], chem_feats[ii])
        ligands.append(lig)

    k_variants = common_k_point_variants(ligands, n_points, min_actives)
    containers = common_k_point_feature_lists(ligands, k_variants)

    pharmacophores = []
    scores = []
    for cont in containers:
        boxes = pharmacophore_partitioning(cont, min_actives)
        for box in boxes:
            box_scores = score_common_pharmacophores(box)
            if len(box_scores) == 0:
                # All pharmacophores exceed RMSD cutoff
                continue
            feat_list = box[box_scores[0][1]]
            lig = ligands[feat_list.id[0]]
            pharma = feat_list_to_pharma(feat_list, lig)

            pharmacophores.append(pharma)
            scores.append(box_scores[0][0])

    return pharmacophores, scores
