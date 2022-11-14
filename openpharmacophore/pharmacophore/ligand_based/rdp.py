import numpy as np
import math

BIN_SIZE = 1.0  # In angstroms
# Maximum interpoint distance in angstroms
MAX_DIST = 15.0
BINS = np.arange(0, MAX_DIST + BIN_SIZE, step=BIN_SIZE)


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
        low_bin = math.floor(flist.distances[dim] - BIN_SIZE)
        upp_bin = low_bin + 1
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
