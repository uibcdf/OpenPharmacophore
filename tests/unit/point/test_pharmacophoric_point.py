from openpharmacophore import PharmacophoricPoint
from openpharmacophore.pharmacophore.exceptions import WrongDimensionalityError, InvalidFeatureError, IncorrectShapeError, \
    NotAQuantityError
import numpy as np
import pyunitwizard as puw
import pytest


@pytest.fixture()
def hydrogen_bond_donor():
    """ Returns a hydrogen bond donor with no direction"""

    radius = puw.quantity(1.0, "angstroms")
    center = puw.quantity([1.0, 1.0, 1.0], "angstroms")

    # Create pharmacophoric point without direction
    atom_inxs = (3, 4, 5, 6)
    return PharmacophoricPoint("hb donor", center, radius, None, atom_inxs)


@pytest.fixture()
def aromatic_ring():
    """ Returns an aromatic ring with no atom indices."""
    # Create pharmacophoric point without atom indices
    center = puw.quantity([1.5, -2.0, 3.2], "angstroms")
    radius = puw.quantity(1.5, "angstroms")
    direction = np.array([1.0, 1.0, 1.0])

    return PharmacophoricPoint("aromatic ring", center, radius, direction)


def test_init_pharmacophoric_point(hydrogen_bond_donor, aromatic_ring):
    # Test hydrogen bond donor
    assert not hydrogen_bond_donor.has_direction
    assert np.allclose(puw.get_value(hydrogen_bond_donor.center, "angstroms"),
                       np.array([1.0, 1.0, 1.0]))
    assert np.allclose(puw.get_value(hydrogen_bond_donor.radius, "angstroms"), 1.0)
    assert hydrogen_bond_donor.atom_indices == {3, 4, 5, 6}
    assert hydrogen_bond_donor.feature_name == "hb donor"

    # Test ring
    ring = aromatic_ring
    assert ring.has_direction
    assert ring.feature_name == "aromatic ring"
    assert len(ring.atom_indices) == 0
    assert np.allclose(puw.get_value(ring.center, "angstroms"),
                       np.array([1.5, -2.0, 3.2]))
    assert np.allclose(puw.get_value(ring.radius, "angstroms"), 1.5)
    assert np.allclose(ring.direction,
                       np.array([[1.0, 1.0, 1.0]]) / np.linalg.norm(np.array([1.0, 1.0, 1.0])))


def test_init_pharmacophoric_point_center_is_not_quantity():
    radius = puw.quantity(1.0, "angstroms")
    with pytest.raises(NotAQuantityError, match="center is of type <class 'list'>"):
        PharmacophoricPoint(feat_type="hb donor", center=[1, 2, 3], radius=radius)


def test_init_pharmacophoric_point_center_has_wrong_dim():
    radius = puw.quantity(1.0, "angstroms")
    with pytest.raises(WrongDimensionalityError, match="center has incorrect dimensionality"):
        PharmacophoricPoint(feat_type="hb donor", center=puw.quantity([1.0, 1.0, 1.0], "seconds"), radius=radius)


def test_init_pharmacophoric_point_center_has_wrong_shape():
    radius = puw.quantity(1.0, "angstroms")
    with pytest.raises(IncorrectShapeError, match="center has incorrect shape"):
        PharmacophoricPoint(feat_type="hb donor", center=puw.quantity([1.0, 1.0], "angstroms"), radius=radius)


def test_init_pharmacophoric_point_radius_is_not_a_quantity():
    center = puw.quantity([1.0, 1.0, 1.0], "angstroms")

    with pytest.raises(NotAQuantityError, match="radius is of type <class 'float'>"):
        PharmacophoricPoint(feat_type="hb donor", center=center, radius=1.0)


def test_init_pharmacophoric_point_radius_has_wrong_dim():
    center = puw.quantity([1.0, 1.0, 1.0], "angstroms")
    with pytest.raises(WrongDimensionalityError, match="radius has incorrect dimensionality"):
        PharmacophoricPoint(feat_type="hb donor", center=center, radius=puw.quantity(1.0, "seconds"))


def test_pharmacophoric_point_init_with_invalid_feature():
    radius = puw.quantity(1.0, "angstroms")
    center = puw.quantity([1.0, 1.0, 1.0], "angstroms")

    with pytest.raises(InvalidFeatureError, match="Invalid feature name rubber duck"):
        PharmacophoricPoint(feat_type="rubber duck", center=center, radius=radius)


def test_pharmacophoric_point_equality_based_on_coordinates():
    radius = puw.quantity(1.0, "angstroms")
    center = puw.quantity([1.0, 1.0, 1.0], "angstroms")
    feat_name = "hb donor"

    donor_1 = PharmacophoricPoint(feat_name, center, radius)
    donor_2 = PharmacophoricPoint(feat_name, center, radius)
    donor_3 = PharmacophoricPoint(feat_name,
                                  puw.quantity([2.0, 2.0, 2.0], "angstroms"),
                                  radius)

    assert donor_1 == donor_2
    assert not donor_1 == donor_3


def test_pharmacophoric_point_string_representation(hydrogen_bond_donor, aromatic_ring):
    expected_donor_str = ("PharmacophoricPoint(feat_type=hb donor; "
                          "center=(1.0, 1.0, 1.0); "
                          "radius=1.0)")

    assert str(hydrogen_bond_donor) == expected_donor_str

    expected_aromatic_str = ("PharmacophoricPoint(feat_type=aromatic ring; "
                             "center=(1.5, -2.0, 3.2); "
                             "radius=1.5; "
                             "direction=(0.58, 0.58, 0.58))")

    assert str(aromatic_ring) == expected_aromatic_str
