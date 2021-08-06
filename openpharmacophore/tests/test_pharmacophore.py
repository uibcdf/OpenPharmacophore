from openpharmacophore.pharmacophore import Pharmacophore
from openpharmacophore import pharmacophoric_elements
from openpharmacophore._private_tools.exceptions import InvalidFeatureError
import pyunitwizard as puw
import pytest

@pytest.fixture
def three_element_pharmacophore():
    """Returns as pharmacophore with three elements"""
    hb_acceptor = pharmacophoric_elements.HBAcceptorSphere(
        center=puw.quantity([1,0,0], "angstroms"),
        radius=puw.quantity(1.0, "angstroms")
    )
    ring_1 = pharmacophoric_elements.AromaticRingSphere(
        center=puw.quantity([2, 1, 4], "angstroms"),
        radius=puw.quantity(1.0, "angstroms")
    )
    ring_2 = pharmacophoric_elements.AromaticRingSphere(
        center=puw.quantity([0, 1, 2], "angstroms"),
        radius=puw.quantity(1.0, "angstroms")
    )
    pharmacophore = Pharmacophore(elements=[hb_acceptor, ring_1, ring_2])
    return pharmacophore

def test_add_element(three_element_pharmacophore):
    hydrophobic = pharmacophoric_elements.HydrophobicSphere(
        center=puw.quantity([-1, 0, 2], "angstroms"),
        radius=puw.quantity(1.0, "angstroms")
    )
    three_element_pharmacophore.add_element(hydrophobic)
    assert three_element_pharmacophore.n_elements == 4
    assert isinstance(three_element_pharmacophore.elements[3], pharmacophoric_elements.HydrophobicSphere)

def test_remove_element(three_element_pharmacophore):
    three_element_pharmacophore.remove_element(0)
    assert three_element_pharmacophore.n_elements == 2

@pytest.mark.parametrize('feat_type,exception', [
    ("aromatic ring", None), ("hydrophobicity", InvalidFeatureError)
])
def test_remove_feature(three_element_pharmacophore, feat_type, exception):
    try:
        three_element_pharmacophore.remove_feature(feat_type)
    except(InvalidFeatureError) as e:
        assert e.message == "Cannot remove feature. The pharmacophore does not contain any hydrophobicity"
    else:
        assert three_element_pharmacophore.n_elements == 1
        assert not isinstance(three_element_pharmacophore.elements[0], pharmacophoric_elements.AromaticRingSphere)

def test_reset(three_element_pharmacophore):
    three_element_pharmacophore._reset()
    assert three_element_pharmacophore.n_elements == 0
    assert len(three_element_pharmacophore.elements) == 0
    assert three_element_pharmacophore.extractor is None
    assert three_element_pharmacophore.molecular_system is None

def test_add_to_NGLView():
    pass

def test_show():
    pass