from openpharmacophore.molecular_systems.ligand import Ligand, LigandWithHsError, smiles_from_pdb_id
import numpy as np
import pyunitwizard as puw
import pytest
from copy import deepcopy


def test_ligand_from_topology(estradiol):
    assert estradiol.n_atoms == 20
    assert estradiol.n_conformers == 1
    assert estradiol.n_bonds == 23
    assert estradiol.lig_id == "EST"


@pytest.fixture
def estradiol_fixed_bonds(estradiol):
    ligand = deepcopy(estradiol)
    ligand.fix_bond_order(
        smiles="C[C@]12CC[C@@H]3c4ccc(cc4CC[C@H]3[C@@H]1CC[C@@H]2O)O"
    )
    return ligand


def test_fix_bond_order(estradiol_fixed_bonds):
    assert estradiol_fixed_bonds.has_aromatic_bonds()
    

@pytest.fixture()
def caffeine():
    return Ligand.from_string(
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        form="smi"
    )


@pytest.fixture()
def caffeine_with_hyd(caffeine):
    ligand = deepcopy(caffeine)
    assert not ligand.has_hydrogens

    ligand.add_hydrogens()
    return ligand


def test_has_hydrogens(caffeine):
    assert not caffeine.has_hydrogens


def test_add_hydrogens_to_ligand_with_no_conformers(caffeine_with_hyd):
    assert caffeine_with_hyd.has_hydrogens
    assert caffeine_with_hyd.n_atoms == 24


def test_add_hydrogens_to_ligand_with_conformers(estradiol_fixed_bonds):
    ligand = estradiol_fixed_bonds
    assert not ligand.has_hydrogens

    ligand.add_hydrogens()
    assert ligand.has_hydrogens
    assert ligand.n_atoms == 44
    assert ligand.get_conformer(0).shape == (44, 3)


def test_add_hydrogens_to_ligand_with_hydrogens_raises_error(caffeine_with_hyd):
    with pytest.raises(LigandWithHsError):
        caffeine_with_hyd.add_hydrogens()


def test_smiles_from_pdb_id_with_chain():
    mapper = {
        "DAO": "CCCCCCCCCCCC(=O)O",
        "EST": "C[C@]12CC[C@@H]3c4ccc(cc4CC[C@H]3[C@@H]1CC[C@@H]2O)O"
    }
    pdb_id = "DAO:B"
    assert smiles_from_pdb_id(pdb_id, mapper=mapper) == "CCCCCCCCCCCC(=O)O"


def test_smiles_from_pdb_id_without_chain():
    mapper = {
        "DAO": "CCCCCCCCCCCC(=O)O",
        "EST": "C[C@]12CC[C@@H]3c4ccc(cc4CC[C@H]3[C@@H]1CC[C@@H]2O)O"
    }
    pdb_id = "DAO"
    assert smiles_from_pdb_id(pdb_id, mapper=mapper) == "CCCCCCCCCCCC(=O)O"


def test_add_conformers(estradiol):
    estradiol = deepcopy(estradiol)
    assert estradiol.n_conformers == 1

    coords = puw.quantity(
        np.ones((2, estradiol.n_atoms, 3)), "angstroms"
    )
    estradiol.add_conformers(coords)
    assert estradiol.n_conformers == 3


def test_get_conformer(estradiol):
    conformer = estradiol.get_conformer(0)
    assert conformer.shape == (20, 3)
    assert str(puw.get_unit(conformer)) == "angstrom"
