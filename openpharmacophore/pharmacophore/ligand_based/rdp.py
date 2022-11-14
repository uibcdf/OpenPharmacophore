import numpy as np

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
    """ A container of FeatureLists"""
    pass


def recursive_partitioning(feat_lists, dim, n_pairs, boxes):
    """ Obtain common pharmacophores by recursive distance partitioning.

        Parameters
        ----------
        feat_lists : FLContainer
        dim : int
        n_pairs : int
        boxes : list[FeatureList]

    """
    bins = [
        FLContainer(bin=(BINS[ii], BINS[ii + 1])) for ii in range(BINS.shape[0] - 1)
    ]
