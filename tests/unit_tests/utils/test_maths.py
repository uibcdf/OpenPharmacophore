import openpharmacophore.utils.maths as maths
import pyunitwizard as puw
import numpy as np


def test_points_distance():
    distance = maths.points_distance(
        puw.quantity(np.array([0., 0., 0.]), "angstroms"),
        puw.quantity(np.array([4., 3., 0.]), "angstroms"),
    )
    assert puw.get_value(distance) == 5.
