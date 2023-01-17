from openpharmacophore.molecular_systems.ligand import Ligand, smiles_from_pdb_id
import numpy as np
import pyunitwizard as puw
import pytest
from copy import deepcopy


def test_ligand_from_topology(estradiol):
    assert estradiol.n_atoms == 20
    assert estradiol.n_conformers == 1
    assert estradiol.n_bonds == 23


def test_fix_bond_order(estradiol):
    ligand = deepcopy(estradiol)
    ligand.fix_bond_order(
        smiles="C[C@]12CC[C@@H]3c4ccc(cc4CC[C@H]3[C@@H]1CC[C@@H]2O)O"
    )
    assert ligand.has_aromatic_bonds()
    

@pytest.fixture()
def caffeine():
    return Ligand.from_string(
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        form="smi"
    )


def test_has_hydrogens(caffeine):
    assert not caffeine.has_hydrogens()


def test_add_hydrogens(caffeine):
    ligand = deepcopy(caffeine)
    ligand.add_hydrogens()
    assert ligand.has_hydrogens()
    assert ligand.n_atoms == 24


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
