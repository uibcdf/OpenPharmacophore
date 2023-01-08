import pytest
import pyunitwizard as puw
from openpharmacophore import PharmacophoricPoint
from pathlib import Path


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
