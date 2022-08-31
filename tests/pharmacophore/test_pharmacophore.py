from openpharmacophore.pharmacophore.pharmacophore import Pharmacophore
from openpharmacophore import LigandBasedPharmacophore
import pytest


def test_cannot_instantiate_pharmacophore():

    with pytest.raises(TypeError):
        ph = Pharmacophore()

# The following tests are called on LigandBasedPharmacophore because we cannot
# instantiate Pharmacophore as is an abstract class


def test_pharmacophore_len():
    assert False, "Implement me!"


def test_iterate_pharmacophore():
    assert False, "Implement me!"


def test_pharmacophore_get_item():
    assert False, "Implement me!"


def test_pharmacophore_equality():
    assert False, "Implement me!"


def test_distance_matrix():
    assert False, "Implement me!"


def test_to_rdkit():
    assert False, "Implement me!"


def test_to_ligand_scout():
    assert False, "Implement me!"


def test_to_moe():
    assert False, "Implement me!"


def test_to_pharmagist():
    assert False, "Implement me!"


def test_pharmacophore_string_representation():
    assert False, "Implement me!"
