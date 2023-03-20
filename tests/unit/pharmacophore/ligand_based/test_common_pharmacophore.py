import pyunitwizard as puw
import numpy as np
from openpharmacophore.pharmacophore.ligand_based.common_pharmacophore import CommonPharmacophoreFinder


def test_init_common_pharmacophore_finder():
    scoring_fn_params = {
        "point_weight": 1.0,
        "rmsd_cutoff": puw.quantity(1.2, "angstroms"),
        "vector_weight": 1.0,
        "cos_cutoff": 0.5
    }
    min_dist = puw.quantity(1, "angstroms")
    max_dist = puw.quantity(10, "angstroms")
    bin_size = min_dist

    finder = CommonPharmacophoreFinder(
        n_points=3,
        min_actives=4,
        max_pharmacophores=10,
        scoring_fn_params=scoring_fn_params,
        max_dist=max_dist,
        min_dist=min_dist,
        bin_size=bin_size,
    )

    expected_bins = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    assert finder.min_dist == min_dist
    assert finder.max_dist == max_dist
    assert finder.bin_size == bin_size
    assert np.all(finder.bins == expected_bins)

    assert finder.scoring_fn.point_weight == 1.0
    assert finder.scoring_fn.vector_weight == 1.0
    assert finder.scoring_fn.rmsd_cutoff == puw.quantity(1.2, "angstroms")
    assert finder.scoring_fn.cos_cutoff == 0.5
