import numpy as np
import math
import itertools

BIN_SIZE = 1.0  # In angstroms
# Maximum interpoint distance in angstroms
MAX_DIST = 15.0
BINS = np.arange(0, MAX_DIST + BIN_SIZE, step=BIN_SIZE)

# Scoring function values
W_POINT = 1.0
W_VECTOR = 1.0
RMSD_CUTOFF = 1.2
COS_CUTOFF = 0.5


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


def common_k_point_variants(variants,  n_points, min_actives):
    """ Find the common k-point variants.

        Parameters
        ----------
        variants : list[str]
            A list with the variant of each ligand. Each variant must include all
            of its ligand chemical features.

        n_points : int


        min_actives : int
            Number of actives that the variant must be present in to be
            considered common

        Returns
        -------
        list[str]
            A list with all the common variants.

    """
    common = {}
    for ii in range(len(variants)):
        for k_var in itertools.combinations(variants[ii], n_points):
            var = "".join(k_var)
            try:
                count = common[var]
                # Only increase count if it is a different ligand
                if count <= ii + 1:
                    common[var] += 1
            except KeyError:
                common[var] = 1

    return [var for var, count in common.items() if count >= min_actives]


def find_common_pharmacophores(ligands, n_points, min_actives):
    """ Find common pharmacophores in a set of ligands and assigns a score to each one.

        Parameters
        ----------
        ligands : list[rdkit.Chem.Mol]
            List of molecules with conformers.

        n_points : int
        min_actives : int

        Returns
        -------

    """
    pass
