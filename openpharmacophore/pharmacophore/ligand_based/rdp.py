from openpharmacophore.pharmacophore.chem_feats import smarts_ligand, feature_indices, feature_centroids
from openpharmacophore.utils.maths import points_distance
import numpy as np
import math
import itertools
from collections import defaultdict, Counter


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

    feature_to_char = {
        "hb acceptor": "A",
        "hb donor": "D",
        "aromatic ring": "R",
        "hydrophobicity": "H",
        "positive charge": "P",
        "negative charge": "N",
    }

    def __init__(self, mol):
        self.mol = mol
        self.feats = {}
        self.variant = ""
        self.feat_count = {}
        self.distances = None

    def find_features(self):
        """ Find chemical features in this ligand.

            Updates the variant, feat_count and distances attributes.

        """
        for feat_type, smarts in smarts_ligand.items():
            indices = feature_indices(smarts, self.mol)
            short_name = self.feature_to_char[feat_type]
            self.feats[short_name] = indices
            self.variant += short_name * len(indices)

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

    def has_k_variant(self, var):
        """ Returns true if the ligand contains a k-variant.

            Parameters
            ----------
            var : str

            Returns
            -------
            bool
        """
        if len(var) <= len(self.variant):
            count = Counter(var)
            for feat, cnt in count.items():
                try:
                    if cnt > self.feat_count[feat]:
                        return False
                except KeyError:
                    return False
            return True
        return False


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

    def __init__(self, variant, fl_id, distances):
        self.n_pairs = len(variant) * (len(variant) - 1) / 2
        if distances.shape[0] != self.n_pairs:
            raise ValueError

        self.id = fl_id
        self.distances = distances
        self.variant = variant
        # TODO: store the indices of the variant in the original ligand


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
        return low_bin, low_bin + 1
    else:
        low_bin = math.ceil(num - bin_size)
        return low_bin, low_bin + 1


def recursive_partitioning(container, dim, n_pairs, boxes, n_mols):
    """ Obtain common pharmacophores by recursive distance partitioning.

        Parameters
        ----------
        container : FLContainer
        dim : int
        n_pairs : int
        boxes : list[FLContainer]
        n_mols : int

    """
    bins = [
        FLContainer(bin=(BINS[ii], BINS[ii + 1])) for ii in range(BINS.shape[0] - 1)
    ]

    for flist in container:
        # Assign each distance to the two closest bins
        # low_bin = math.floor(flist.distances[dim] - BIN_SIZE)
        # upp_bin = low_bin + 1
        low_bin, upp_bin = nearest_bins(flist.distances[dim], BIN_SIZE)
        bins[low_bin].append(flist)
        bins[upp_bin].append(flist)

    for container in bins:
        # We process a bin only if it contains at least one conformer
        # from each molecule
        if len(container.mols) == n_mols:
            if dim < n_pairs - 1:
                recursive_partitioning(
                    container, dim + 1, n_pairs, boxes, n_mols
                )
            else:
                boxes.append(container)


def score_common_pharmacophores(box):
    """ Calculate the score of all feature lists in a container and sort them by score.

        Parameters
        ----------
        box : FLContainer

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
            rmsd = np.sqrt(((box[ii].distances - box[jj].distances) ** 2).mean())
            if rmsd >= RMSD_CUTOFF:
                exclude.append(ii)
            point_score = 1 - rmsd / RMSD_CUTOFF
            vec_score = vector_score(box[ii].id, box[jj].id)
            score = W_POINT * point_score + W_VECTOR * vec_score
            scores.append((score, ii, jj))

    scores.sort(reverse=True)
    if len(exclude) > 0:
        # Filter pharmacophores that exceed rmsd cutoff
        return [s for s in scores if s[1] not in exclude]
    else:
        return scores


def vector_score(id_1, id_2):
    # TODO: complete me!
    return 0


def common_k_point_variants(ligands,  n_points, min_actives):
    """ Find the common k-point variants, that is, the variants consisting
        of k pharmacophoric points that are common to at least the specified
        number of actives.

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
        list[str]
            A list with all the common k-point variants.

    """
    common = defaultdict(int)
    for ii in range(len(ligands)):
        # Keep track of the variants in this ligand
        mol_variant = {}
        for k_var in itertools.combinations(ligands[ii].variant, n_points):
            var = "".join(k_var)
            try:
                mol_variant[var] += 1
            except KeyError:
                common[var] += 1
                mol_variant[var] = 1

    return [var for var, count in common.items() if count >= min_actives]


def common_k_point_feature_lists(ligands, k_variants):
    """ Returns a feature list container for each of the k-point variants.

        Parameters
        ----------
        ligands : list[Ligand]
        k_variants : list[str]

        Returns
        -------
        list[FLContainer]

    """
    all_containers = {
        var: FLContainer(var) for var in k_variants
    }
    for var in k_variants:
        for ii, lig in enumerate(ligands):
            if lig.has_k_variant(var):
                # TODO: if a variant has repeated feature type. We must create a feat_list
                #  for each possible combination of this repeated features.
                container = all_containers[var]
                for jj in range(lig.num_confs):
                    fl_id = (ii, jj)
                    distances = lig.k_distances(jj, var)
                    feat_list = FeatureList(var, fl_id, distances)
                    container.append(feat_list)

    return list(all_containers.values())


def find_common_pharmacophores(mols, n_points, min_actives):
    """ Find common pharmacophores in a set of ligands and assigns a score to each one.

        Parameters
        ----------
        mols : list[rdkit.Chem.Mol]
            List of molecules with conformers.

        n_points : int
            Number of pharmacophoric points the common pharmacophores will have.

        min_actives : int
            Number of ligands that the common pharmacophores are present in.

        Returns
        -------

    """
    ligands = []
    for mol in mols:
        lig = Ligand(mol)
        lig.find_features()
        ligands.append(lig)

    k_variants = common_k_point_variants(ligands, n_points, min_actives)
    containers = common_k_point_feature_lists(ligands, k_variants)

    boxes = []
    for cont in containers:
        box = []
        recursive_partitioning(cont, 0, cont[0].n_pairs, box, min_actives)
        boxes.append(box)

    cps = []
    for box in boxes:
        cps.append(score_common_pharmacophores(box)[0])
