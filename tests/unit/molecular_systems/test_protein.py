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
