from openpharmacophore import PharmacophoricPoint
from openpharmacophore.pharmacophore.pharmacophoric_point import distance_bewteen_pharmacophoric_points
import openpharmacophore._private_tools.exceptions as exc
import numpy as np
import pyunitwizard as puw
import pytest


def test_pharmacophoric_point_constructor():

    # First we test that it raises the correct exceptions when the input arguments
    # are not valid

    radius = puw.quantity(1.0, "angstroms")
    center = puw.quantity([1.0, 1.0, 1.0], "angstroms")
    
    # Test that center is validated correctly
    with pytest.raises(exc.IsNotQuantityError, match="center is not a quantity"):
        PharmacophoricPoint(feat_type="hb donor", center=[1, 2, 3], radius=radius)
    with pytest.raises(exc.WrongDimensionalityError, match="center has dimensionality"):
        PharmacophoricPoint(feat_type="hb donor", center=puw.quantity([1.0, 1.0, 1.0], "seconds"), radius=radius)    
    with pytest.raises(exc.BadShapeError, match="center has shape"):
        PharmacophoricPoint(feat_type="hb donor", center=puw.quantity([1.0, 1.0], "angstroms"), radius=radius)
    # Test for radius
    with pytest.raises(exc.IsNotQuantityError, match="radius is not a quantity"):
        PharmacophoricPoint(feat_type="hb donor", center=center, radius=1.0)
    with pytest.raises(exc.NegativeRadiusError, match="radius must be a positive"):
        PharmacophoricPoint(feat_type="hb donor", center=center, radius=puw.quantity(-1.0, "angstroms"))
    with pytest.raises(exc.WrongDimensionalityError, match="radius has dimensionality"):
        PharmacophoricPoint(feat_type="hb donor", center=center, radius=puw.quantity(1.0, "seconds")) 
    # Test feature type
    with pytest.raises(exc.InvalidFeatureType, match="is not a valid feature type"):
        PharmacophoricPoint(feat_type="rubber duck", center=center, radius=radius)
    with pytest.raises(exc.OpenPharmacophoreTypeError, match="atom_indices must be"):
        PharmacophoricPoint(feat_type="hb donor", center=center, radius=radius, atom_indices=1)

    feat_name = "hb donor"
    atom_inxs = (3, 4, 5, 6)

    donor_1 = PharmacophoricPoint(feat_name, center, radius, None, atom_inxs)
    donor_2 = PharmacophoricPoint(feat_name, center, radius, 
                direction=np.array([1.0, 1.0, 1.0]), atom_indices=atom_inxs)

    assert donor_1.has_direction == False
    assert np.allclose(puw.get_value(donor_1.center, "angstroms"), np.array([1.0, 1.0, 1.0]))
    assert np.allclose(puw.get_value(donor_1.radius, "angstroms"), 1.0)
    assert atom_inxs == (3, 4, 5, 6) 

    feat_name = "aromatic ring"
    center = puw.quantity([1.5, -2.0, 3.2], "angstroms")
    radius = puw.quantity(1.5, "angstroms")
    direction = np.array([1.0, 1.0, 1.0])
    atom_inxs = None

    ring = PharmacophoricPoint(feat_name, center, radius, direction, atom_inxs)
    assert ring.has_direction == True
    assert ring.atom_indices is None
    assert np.allclose(puw.get_value(ring.center, "angstroms"), np.array([1.5, -2.0, 3.2]))
    assert np.allclose(puw.get_value(ring.radius, "angstroms"), 1.5) 
    assert np.allclose(ring.direction, 
            np.array([[1.0, 1.0, 1.0]]) / np.linalg.norm(np.array([1.0, 1.0, 1.0])))

    assert donor_1 != ring
    assert donor_2 != donor_1


def test_distance_between_pharmacophoric_points():
    
    radius = puw.quantity(1.0, "angstroms") 

    point_1 = PharmacophoricPoint("aromatic ring", puw.quantity([3, 4, 0], "angstroms"), radius)
    point_2 = PharmacophoricPoint("hb acceptor", puw.quantity([0, 0, 0], "angstroms"), radius)
    
    assert distance_bewteen_pharmacophoric_points(point_1, point_2) == 5

    point_1 = PharmacophoricPoint("aromatic ring", puw.quantity([1, 1, 1], "angstroms"), radius)
    point_2 = PharmacophoricPoint("hb acceptor", puw.quantity([1, 1, 1], "angstroms"), radius)
    
    assert distance_bewteen_pharmacophoric_points(point_1, point_2) == 0


def test_add_to_ngl_view():
    assert False, "Complete me!"


def test_pharmacophoric_point_string_representation():
    pass
