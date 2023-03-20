from openpharmacophore import PharmacophoricPoint, distance_between_pharmacophoric_points
from openpharmacophore.point.exceptions import WrongDimensionalityError, InvalidFeatureError, IncorrectShapeError, \
    NotAQuantityError
import numpy as np
import nglview as nv
import pyunitwizard as puw
import pytest
from matplotlib.colors import to_rgb
from copy import deepcopy
from unittest.mock import Mock, call


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


def test_add_to_ngl_view_point_with_no_direction(hydrogen_bond_donor):

    mock_view = Mock()
    mock_view._ngl_component_ids = []
    hydrogen_bond_donor.add_to_ngl_view(mock_view)
    donor_color = to_rgb("#17A589")
    radius = 1.0
    center = puw.get_value(hydrogen_bond_donor.center, "angstroms").tolist()

    assert mock_view.shape.add_sphere.call_args == call(
        center, donor_color, radius, "hb donor"
    )
    assert mock_view.update_representation.call_args == call(
        component=0, repr_index=0, opacity=0.5
    )


def test_add_to_ngl_view_point_with_direction(aromatic_ring):
    mock_view = Mock()
    mock_view._ngl_component_ids = [0, 1]
    aromatic_ring.add_to_ngl_view(mock_view, opacity=0.7)

    ring_color = to_rgb("#F1C40F")
    radius = 1.5
    center = puw.get_value(aromatic_ring.center, "angstroms").tolist()
    arrow_end = (center + radius * aromatic_ring.direction*2.0).tolist()

    assert mock_view.shape.add_sphere.call_args == call(
        center, ring_color, radius, "aromatic ring"
    )
    assert mock_view.shape.add_arrow.call_args == call(
        center, arrow_end, ring_color, 0.2
    )
    assert mock_view.update_representation.call_count == 2
    assert mock_view.update_representation.call_args == call(
        component=3, repr_index=0, opacity=0.9
    )


def test_add_to_ngl_view_components_are_loaded(hydrogen_bond_donor, aromatic_ring):

    view = nv.NGLWidget()
    assert len(view._ngl_component_ids) == 0

    # HB donor with no direction should just add a sphere to the visualization
    hydrogen_bond_donor.add_to_ngl_view(view)
    assert len(view._ngl_component_ids) == 1
    assert view._ngl_component_names[0] == "nglview.shape.Shape"

    # Aromatic ring should add a sphere and an arrow to the visualization
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


def test_set_center_pharmacophoric_point(hydrogen_bond_donor):
    donor = deepcopy(hydrogen_bond_donor)
    donor.center = puw.quantity(np.array([3.5, 2.2, 1.0]), "angstroms")
    assert np.all(puw.get_value(donor.center) == np.array([3.5, 2.2, 1.0]))


def test_distance_between_pharmacophoric_points():
    radius = puw.quantity(1.0, "angstroms")

    point_1 = PharmacophoricPoint("aromatic ring", puw.quantity([3, 4, 0], "angstroms"), radius)
    point_2 = PharmacophoricPoint("hb acceptor", puw.quantity([0, 0, 0], "angstroms"), radius)

    assert distance_between_pharmacophoric_points(point_1, point_2) == 5

    point_1 = PharmacophoricPoint("aromatic ring", puw.quantity([1, 1, 1], "angstroms"), radius)
    point_2 = PharmacophoricPoint("hb acceptor", puw.quantity([1, 1, 1], "angstroms"), radius)

    assert distance_between_pharmacophoric_points(point_1, point_2) == 0
