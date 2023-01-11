import pyunitwizard as puw
import numpy as np
from copy import deepcopy

from openpharmacophore import Protein


def test_remove_ligand(topology_with_ligand):
    topology = deepcopy(topology_with_ligand)
    coords = puw.quantity(
        np.ones((2, topology.n_atoms, 3)),
        "nanometers"
    )
    protein = Protein(topology, coords)
    n_atoms = protein.n_atoms

    protein._remove_ligand_by_indices([4, 5, 6, 7])
    assert protein.n_atoms == n_atoms - 4
    assert protein._coords.shape == (2, n_atoms - 4, 3)

    expected_coords = np.ones((2, n_atoms - 4, 3))
    assert np.all(puw.get_value(protein._coords, "nanometers")
                  == expected_coords)


def test_residues_at_distance(topology_2_chains):
    centroid = puw.quantity(np.zeros((1, 3)), "nanometers")
    max_dist = puw.quantity(3.0, "nanometers")

    coords = puw.quantity(np.array([[
        [4., 4., 4.],
        [1., 1., 1.],
        [4., 4., 4.],
        [1., 1., 1.],
        [4., 4., 4.],
        [4., 4., 4.],
        [4., 4., 4.],
        [4., 4., 4.],
        [1., 1., 1.],
    ]]), "nanometers")
    protein = Protein(topology_2_chains, coords)
    indices = protein.residues_at_distance(0, centroid, max_dist)
    expected = np.array([1, 3, 8])
    assert np.all(indices == expected)
