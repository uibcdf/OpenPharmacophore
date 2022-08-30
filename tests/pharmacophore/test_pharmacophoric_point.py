from openpharmacophore import PharmacophoricPoint, distance_between_pharmacophoric_points
import openpharmacophore._private_tools.exceptions as exc
import numpy as np
import nglview as nv
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
    assert not hydrogen_bond_donor.has_direction
    assert np.allclose(puw.get_value(hydrogen_bond_donor.center, "angstroms"),
                       np.array([1.0, 1.0, 1.0]))
    assert np.allclose(puw.get_value(hydrogen_bond_donor.radius, "angstroms"), 1.0)
    assert hydrogen_bond_donor.atom_indices == {3, 4, 5, 6}
    assert hydrogen_bond_donor.feature_name == "hb donor"

    ring = aromatic_ring
    assert ring.has_direction
    assert ring.feature_name == "aromatic ring"
    assert len(ring.atom_indices) == 0
    assert np.allclose(puw.get_value(ring.center, "angstroms"),
                       np.array([1.5, -2.0, 3.2]))
    assert np.allclose(puw.get_value(ring.radius, "angstroms"), 1.5)
    assert np.allclose(ring.direction,
                       np.array([[1.0, 1.0, 1.0]]) / np.linalg.norm(np.array([1.0, 1.0, 1.0])))


def test_pharmacophoric_center_validation():
    # First we test that it raises the correct exceptions when the input arguments
    # are not valid
    radius = puw.quantity(1.0, "angstroms")

    with pytest.raises(exc.IsNotQuantityError, match="center is not a quantity"):
        PharmacophoricPoint(feat_type="hb donor", center=[1, 2, 3], radius=radius)
    with pytest.raises(exc.WrongDimensionalityError, match="center has dimensionality"):
        PharmacophoricPoint(feat_type="hb donor", center=puw.quantity([1.0, 1.0, 1.0], "seconds"), radius=radius)
    with pytest.raises(exc.BadShapeError, match="center has shape"):
        PharmacophoricPoint(feat_type="hb donor", center=puw.quantity([1.0, 1.0], "angstroms"), radius=radius)


def test_pharmacophoric_point_radius_validation():
    center = puw.quantity([1.0, 1.0, 1.0], "angstroms")

    with pytest.raises(exc.IsNotQuantityError, match="radius is not a quantity"):
        PharmacophoricPoint(feat_type="hb donor", center=center, radius=1.0)
    with pytest.raises(exc.NegativeRadiusError, match="radius must be a positive"):
        PharmacophoricPoint(feat_type="hb donor", center=center, radius=puw.quantity(-1.0, "angstroms"))
    with pytest.raises(exc.WrongDimensionalityError, match="radius has dimensionality"):
        PharmacophoricPoint(feat_type="hb donor", center=center, radius=puw.quantity(1.0, "seconds"))


def test_pharmacophoric_point_init_with_invalid_feature():
    radius = puw.quantity(1.0, "angstroms")
    center = puw.quantity([1.0, 1.0, 1.0], "angstroms")

    with pytest.raises(exc.InvalidFeatureType, match="is not a valid feature type"):
        PharmacophoricPoint(feat_type="rubber duck", center=center, radius=radius)


def test_pharmacophoric_point_atom_indices_validation():
    radius = puw.quantity(1.0, "angstroms")
    center = puw.quantity([1.0, 1.0, 1.0], "angstroms")

    with pytest.raises(exc.NotArrayLikeError, match="atom_indices must be"):
        PharmacophoricPoint(feat_type="hb donor", center=center, radius=radius, atom_indices=1)


def test_pharmacophoric_point_equality_based_on_atom_indices(hydrogen_bond_donor):
    radius = puw.quantity(1.0, "angstroms")
    center = puw.quantity([1.0, 1.0, 1.0], "angstroms")
    feat_name = "hb donor"

    donor_1 = hydrogen_bond_donor
    donor_2 = PharmacophoricPoint(feat_name, center, radius,
                                  direction=np.array([1.0, 1.0, 1.0]),
                                  atom_indices=[2, 3, 4])

    assert not donor_1.is_equal(donor_2)

    donor_3 = PharmacophoricPoint(feat_name, center, radius,
                                  direction=np.array([2.0, 1.0, 1.0]),
                                  atom_indices=[3, 4, 5, 6])

    assert donor_1.is_equal(donor_3)


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


def test_distance_between_pharmacophoric_points():
    radius = puw.quantity(1.0, "angstroms")

    point_1 = PharmacophoricPoint("aromatic ring", puw.quantity([3, 4, 0], "angstroms"), radius)
    point_2 = PharmacophoricPoint("hb acceptor", puw.quantity([0, 0, 0], "angstroms"), radius)

    assert distance_between_pharmacophoric_points(point_1, point_2) == 5

    point_1 = PharmacophoricPoint("aromatic ring", puw.quantity([1, 1, 1], "angstroms"), radius)
    point_2 = PharmacophoricPoint("hb acceptor", puw.quantity([1, 1, 1], "angstroms"), radius)

    assert distance_between_pharmacophoric_points(point_1, point_2) == 0


def test_add_to_ngl_view(hydrogen_bond_donor, aromatic_ring):

    view = nv.NGLWidget()
    assert len(view._ngl_component_ids) == 0

    # HB donor with no direction should just add a sphere to the view
    hydrogen_bond_donor.add_to_ngl_view(view)
    assert len(view._ngl_component_ids) == 1
    assert view._ngl_component_names[0] == "nglview.shape.Shape"

    # Aromatic ring should add a sphere and an arrow to the view
    aromatic_ring.add_to_ngl_view(view)
    assert len(view._ngl_component_ids) == 3
    assert view._ngl_component_names[1] == "nglview.shape.Shape"
    assert view._ngl_component_names[2] == "nglview.shape.Shape"


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
