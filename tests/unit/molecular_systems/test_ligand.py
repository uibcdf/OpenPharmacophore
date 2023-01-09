from openpharmacophore.molecular_systems import Topology
from openpharmacophore.molecular_systems.ligand import Ligand, ligand_from_topology
import numpy as np
import pyunitwizard as puw
import pytest


@pytest.fixture
def estradiol_topology():
    top = Topology()
    top.set_num_chains(1)
    top.add_atoms_to_chain({
        "EST": [
            ("C1", "C"),
            ("C2", "C"),
            ("C3", "C"),
            ("O3", "O"),
            ("C4", "C"),
            ("C5", "C"),
            ("C6", "C"),
            ("C7", "C"),
            ("C8", "C"),
            ("C9", "C"),
            ("C10", "C"),
            ("C11", "C"),
            ("C12", "C"),
            ("C13", "C"),
            ("C14", "C"),
            ("C15", "C"),
            ("C16", "C"),
            ("C17", "C"),
            ("O17", "O"),
            ("C18", "C"),
        ]
    }, 0)
    top.add_bonds_from_dict({
        0: [1, 10],
        1: [0, 2],
        2: [1, 3, 4],
        3: [2, ],
        4: [2, 5, ],
        5: [4, 6, 10],
        6: [5, 7],
        7: [6, 8],
        8: [7, 9, 14],
        9: [8, 10, 11],
        10: [0, 5, 9],
        11: [9, 12],
        12: [11, 13],
        13: [12, 14, 17, 19],
        14: [8, 13, 15],
        15: [14, 16],
        16: [15, 17],
        17: [13, 16, 18],
        18: [17],
        19: [13],
    })
    return top


@pytest.fixture
def estradiol_coords():
    return puw.quantity(np.array([
        [104.106,  17.203,  24.775],
        [102.995,  17.834,  25.370],
        [101.695,  17.355,  25.120],
        [100.598,  17.990,  25.704],
        [101.506,  16.240,  24.274],
        [102.621,  15.588,  23.660],
        [102.371,  14.379,  22.735],
        [103.644,  13.753,  22.086],
        [104.898,  13.873,  22.953],
        [105.178,  15.388,  23.261],
        [103.957,  16.078,  23.918],
        [106.462,  15.459,  24.125],
        [107.711,  14.803,  23.508],
        [107.463,  13.343,  23.124],
        [106.170,  13.270,  22.242],
        [106.228,  11.821,  21.792],
        [107.701,  11.713,  21.263],
        [108.494,  12.719,  22.135],
        [109.610,  12.027,  22.746],
        [107.379,  12.449,  24.419],
    ]), "angstroms")


def test_ligand_from_topology(estradiol_topology,
                              estradiol_coords):
    ligand = ligand_from_topology(estradiol_topology, estradiol_coords)
    assert ligand.n_atoms == 20
    assert ligand.n_conformers == 1
    assert ligand.n_bonds == 23

