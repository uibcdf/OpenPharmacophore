from openpharmacophore.molecular_systems import Topology
from openpharmacophore.molecular_systems.topology_to_mol import topology_to_mol
import pytest
import numpy as np


@pytest.fixture
def dipeptide_topology():
    topology = Topology()
    topology.set_num_chains(1)
    topology.add_atoms_to_chain({
        "THR": [
            ("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O"),
            ("CB", "C"), ("OG1", "O"), ("CG2", "C"),
        ],
        "TYR": [
            ("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O"),
            ("CB", "C"), ("CG", "C"), ("CD1", "C"), ("CD2", "C"),
            ("CE1", "C"), ("CE2", "C"), ("CZ", "C"), ("OH", "O"),
        ],
    }, 0)
    return topology


@pytest.fixture
def dipeptide_coords():
    return np.array([[
        [4.4235,  8.0308,  1.8419],
        [4.3549,  7.9243,  1.7706],
        [4.4528,  7.8252,  1.7077],
        [4.5699,  7.8559,  1.6853],
        [4.2608,  7.9792,  1.6611],
        [4.3375,  8.0468,  1.5608],
        [4.1586,  8.0762,  1.7208],
        [4.4030,  7.7052,  1.6814],
        [4.4799,  7.6021,  1.6156],
        [4.4189,  7.5782,  1.4791],
        [4.3007,  7.6033,  1.4583],
        [4.4721,  7.4725,  1.6954],
        [4.5253,  7.4836,  1.8355],
        [4.6598,  7.4629,  1.8621],
        [4.4412,  7.5147,  1.9416],
        [4.7098,  7.4729,  1.9906],
        [4.4895,  7.5245,  2.0707],
        [4.6246,  7.5036,  2.0946],
        [4.6748,  7.5126,  2.2224],
    ]])


def test_topology_to_mol(dipeptide_topology, dipeptide_coords):
    mol = topology_to_mol(
        dipeptide_topology, dipeptide_coords, remove_hyd=True
    )
    assert mol.GetNumAtoms() == 19
