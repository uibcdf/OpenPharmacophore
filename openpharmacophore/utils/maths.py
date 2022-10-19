import numpy as np


def points_distance(coords_1, coords_2):
    """ Returns the distance between two points in 3D space.

        Parameters
        ----------
        coords_1 : puw.Quantity
        coords_2 : puw.Quantity

        Returns
        -------
        puw.Quantity
    """
    return np.sqrt(np.sum(np.power(coords_1 - coords_2, 2)))
