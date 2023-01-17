import pytest
import pyunitwizard as puw
import numpy as np
from pathlib import Path

from openpharmacophore import PharmacophoricPoint
from openpharmacophore.molecular_systems import Topology
from openpharmacophore.molecular_systems.ligand import ligand_from_topology


# TODO: our unit tests will run faster if we do not use files
data_path = Path(__file__).parent / "data"
pharmacophores_path = data_path / "pharmacophores"
ligands_path = data_path / "ligands"
pdb_path = data_path / "pdb"
traj_path = data_path / "traj"


@pytest.fixture
def json_pharmacophore_path():
    return str(pharmacophores_path / "pharmer/1M70.json")


@pytest.fixture
def ligand_scout_pharmacophore_path():
    return str(pharmacophores_path / "ligandscout/ligscout.pml")


@pytest.fixture
def moe_pharmacophore_path():
    return str(pharmacophores_path / "moe/gmp.ph4")


@pytest.fixture
def mol2_pharmacophore_path_elastase():
    return str(pharmacophores_path / "mol2/elastase.mol2")


@pytest.fixture
def mol2_pharmacophore_path_streptadivin():
    return str(pharmacophores_path / "mol2/streptadivin.mol2")


@pytest.fixture
def ligands_mol2():
    return str(ligands_path / "ace.mol2")


@pytest.fixture()
def ligands_smi():
    return str(ligands_path / "clique_detection.smi")


@pytest.fixture
def ligands_sdf():
    return str(ligands_path / "sdf_example.sdf")


@pytest.fixture
def pdb_1ncr_path():
    return str(pdb_path / "1ncr.pdb")


@pytest.fixture
def small_pdb_with_ligand():
    return str(pdb_path / "test_with_lig.pdb")


@pytest.fixture
def small_pdb_with_no_ligand_1():
    return str(pdb_path / "test_no_lig.pdb")


@pytest.fixture
def small_pdb_with_no_ligand_2():
    return str(pdb_path / "test_no_lig_2.pdb")


@pytest.fixture
def estradiol_pdb():
    return str(pdb_path / "estradiol.pdb")


@pytest.fixture
def small_trajectory_path():
    return str(traj_path / "pentalanine_small.gro")


@pytest.fixture
def ligand_trajectory():
    return str(traj_path / "ligand_traj.gro")


@pytest.fixture
def two_element_pharmacophore():
    """Returns a pharmacophore with an aromatic ring and a hb acceptor"""
    radius = puw.quantity(1.0, "angstroms")
    ring = PharmacophoricPoint(
        feat_type="aromatic ring",
        center=puw.quantity([1, 0, 0], "angstroms"),
        radius=radius,
        direction=[0, 0, 1],
    )
    acceptor = PharmacophoricPoint(
        feat_type="hb acceptor",
        center=puw.quantity([1, 2, 2], "angstroms"),
        direction=[0, 1, 1],
        radius=radius)
    return [ring, acceptor]


@pytest.fixture
def three_element_pharmacophore():
    radius = puw.quantity(1.0, "angstroms")
    ring = PharmacophoricPoint(
        feat_type="aromatic ring",
        center=puw.quantity([1, 0, 0], "angstroms"),
        radius=radius,
        direction=[0, 0, 1],
    )
    acceptor = PharmacophoricPoint(
        feat_type="hb acceptor",
        center=puw.quantity([1, 2, 2], "angstroms"),
        radius=radius,
        direction=[0, 1, 1]
    )
    excluded = PharmacophoricPoint(
        feat_type="excluded volume",
        center=puw.quantity([2, 1, 2], "angstroms"),
        radius=radius)
    return [acceptor, excluded, ring]


@pytest.fixture
def five_element_pharmacophore():
    radius = puw.quantity(1.0, "angstroms")
    ring_1 = PharmacophoricPoint(
        feat_type="aromatic ring",
        center=puw.quantity([1, 0, 0], "angstroms"),
        radius=radius,
        direction=[0, 0, 1],
    )
    ring_2 = PharmacophoricPoint(
        feat_type="aromatic ring",
        center=puw.quantity([0, 1, 2], "angstroms"),
        radius=puw.quantity(1.0, "angstroms"),
        direction=[0, 0, 1],
    )
    acceptor = PharmacophoricPoint(
        feat_type="hb acceptor",
        center=puw.quantity([1, 2, 2], "angstroms"),
        radius=radius,
        direction=[0, 1, 1],
    )
    hb_donor = PharmacophoricPoint(
        feat_type="hb donor",
        center=puw.quantity([1, 2, 2], "angstroms"),
        radius=radius,
        direction=[0, 1, 1]
    )
    hydrophobicity = PharmacophoricPoint(
        feat_type="hydrophobicity",
        center=puw.quantity([-1, 2, 2], "angstroms"),
        radius=radius,
    )
    return [acceptor, hb_donor,
            hydrophobicity, ring_1, ring_2]


@pytest.fixture
def topology_2_chains():
    top = Topology()
    top.set_num_chains(2)
    top.add_atoms_to_chain(
        {
            "ALA": [("CA", "C"), ("O", "O"), ("NA", "N")],
            "MET": [("CA", "C"), ("O", "O"), ("CB", "C")],
        },
        0
    )
    top.add_atoms_to_chain(
        {
            "PRO": [("CA", "C"), ("O", "O"), ("NA", "N")],
            "LEU": [("CA", "C"), ("O", "O"), ("CB", "C")],
        },
        1
    )
    assert top.n_residues == 4
    assert top.n_atoms == 12
    return top


@pytest.fixture()
def topology_with_ligand():
    top = Topology()
    top.set_num_chains(5)
    top.add_atoms_to_chain(
        {
            "ALA": [("CA", "C"), ("O", "O"), ("NA", "N"), ("H", "H")],
        },
        0
    )
    top.add_atoms_to_chain(
        {
            "EST": [("C", "C"), ("O", "O"), ("C", "C"), ("H", "H")],
        },
        1
    )
    top.add_atoms_to_chain(
        {
            "HOH": [("O", "O"), ("H", "H"), ("H", "H")],
        },
        2
    )
    top.add_atoms_to_chain(
        {
            "SO4": [("S", "S"), ("O", "O"), ("O", "O"), ("O", "O"), ("O", "O")],
        },
        3
    )
    top.add_atoms_to_chain(
        {
            "EST": [("C", "C"), ("O", "O"), ("C", "C"), ("H", "H")],
        },
        4
    )
    return top


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
    coords = puw.quantity(np.array([[
        [104.106, 17.203, 24.775],
        [102.995, 17.834, 25.370],
        [101.695, 17.355, 25.120],
        [100.598, 17.990, 25.704],
        [101.506, 16.240, 24.274],
        [102.621, 15.588, 23.660],
        [102.371, 14.379, 22.735],
        [103.644, 13.753, 22.086],
        [104.898, 13.873, 22.953],
        [105.178, 15.388, 23.261],
        [103.957, 16.078, 23.918],
        [106.462, 15.459, 24.125],
        [107.711, 14.803, 23.508],
        [107.463, 13.343, 23.124],
        [106.170, 13.270, 22.242],
        [106.228, 11.821, 21.792],
        [107.701, 11.713, 21.263],
        [108.494, 12.719, 22.135],
        [109.610, 12.027, 22.746],
        [107.379, 12.449, 24.419],
    ]]), "angstroms")
    assert coords.shape == (1, 20, 3)
    return coords


@pytest.fixture
def estradiol(estradiol_topology, estradiol_coords):
    return ligand_from_topology(estradiol_topology, estradiol_coords)


# Data for integration tests

@pytest.fixture
def thrombin_ligands():
    return str(ligands_path / "thrombin_ligands.sdf")


@pytest.fixture
def pdb_3bbh_with_hydrogen():
    return str(pdb_path / "3bbh_hyd.pdb")


@pytest.fixture
def pdb_1m7w():
    return str(pdb_path / "1m7w_A_chain.pdb")


@pytest.fixture
def pdb_er_alpha():
    return str(pdb_path / "er_alpha_A_chain.pdb")


@pytest.fixture
def traj_er_alpha():
    return str(traj_path / "eralpha_small.h5")
