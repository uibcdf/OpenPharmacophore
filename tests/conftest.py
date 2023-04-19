from rdkit import Chem
import pytest
import pyunitwizard as puw
import numpy as np
from pathlib import Path

from openpharmacophore import PharmacophoricPoint, Pharmacophore
from openpharmacophore.molecular_systems import Topology, Protein
from openpharmacophore.molecular_systems.ligand import ligand_from_topology


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
def mol2_pharmacophore_path():
    return str(pharmacophores_path / "mol2/pharmacophore.mol2")


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
    return Pharmacophore([ring, acceptor])


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
    return Pharmacophore([acceptor, excluded, ring])


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
    return Pharmacophore([
        acceptor, hb_donor,
        hydrophobicity, ring_1, ring_2]
    )


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


@pytest.fixture
def protein_4_residues(topology_2_chains):
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
        [4., 4., 4.],
        [4., 4., 4.],
        [4., 4., 4.],
    ]]), "nanometers")
    return Protein(topology_2_chains, coords)


@pytest.fixture
def dialanine():
    # Returns two alanine residues in a single molecule
    mol = Chem.MolFromPDBBlock(
        """MODEL    0
ATOM      1  N   ALA A   1      46.411  13.580  27.508  1.00  0.00           N  
ATOM      2  CA  ALA A   1      46.204  13.395  26.090  1.00  0.00           C  
ATOM      3  C   ALA A   1      47.270  14.190  25.332  1.00  0.00           C  
ATOM      4  O   ALA A   1      47.947  13.654  24.448  1.00  0.00           O  
ATOM      5  CB  ALA A   1      44.817  13.873  25.699  1.00  0.00           C  
ATOM      6  N   ALA A   2      50.424  13.298  22.736  1.00  0.00           N  
ATOM      7  CA  ALA A   2      50.417  13.697  21.333  1.00  0.00           C  
ATOM      8  C   ALA A   2      51.805  14.263  21.035  1.00  0.00           C  
ATOM      9  O   ALA A   2      52.408  13.939  20.017  1.00  0.00           O  
ATOM     10  CB  ALA A   2      49.360  14.751  21.095  1.00  0.00           C  
ATOM     11  H1  ALA A   1      45.750  12.989  28.032  1.00  0.00           H  
ATOM     12  H2  ALA A   1      47.376  13.314  27.752  1.00  0.00           H  
ATOM     13  H3  ALA A   1      46.286  12.328  25.836  1.00  0.00           H  
ATOM     14  H4  ALA A   1      47.435  15.250  25.577  1.00  0.00           H  
ATOM     15  H5  ALA A   1      44.230  14.075  26.607  1.00  0.00           H  
ATOM     16  H6  ALA A   1      44.901  14.794  25.104  1.00  0.00           H  
ATOM     17  H7  ALA A   1      44.316  13.096  25.104  1.00  0.00           H  
ATOM     18  H1  ALA A   2      50.409  14.138  23.332  1.00  0.00           H  
ATOM     19  H2  ALA A   2      49.593  12.722  22.933  1.00  0.00           H  
ATOM     20  H3  ALA A   2      50.187  12.844  20.678  1.00  0.00           H  
ATOM     21  H4  ALA A   2      52.268  14.968  21.741  1.00  0.00           H  
ATOM     22  H5  ALA A   2      48.988  15.122  22.061  1.00  0.00           H  
ATOM     23  H6  ALA A   2      49.796  15.585  20.525  1.00  0.00           H  
ATOM     24  H7  ALA A   2      48.527  14.313  20.525  1.00  0.00           H  
END
""", removeHs=False)
    assert mol.GetNumAtoms() == 24
    return mol


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
