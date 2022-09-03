from openpharmacophore import LigandBasedPharmacophore, PharmacophoricPoint
import numpy as np
import pyunitwizard as puw
import pytest


def test_init_ligand_based_pharmacophore():
    lbp = LigandBasedPharmacophore(None)


@pytest.fixture()
def pharmacophore_three_points():
    radius = puw.quantity(1.0, "angstroms")
    center = puw.quantity([1.0, 1.0, 1.0], "angstroms")
    points = [
        PharmacophoricPoint("hb donor", center, radius),
        PharmacophoricPoint("aromatic ring", center * 2.0, radius),
        PharmacophoricPoint("hydrophobicity", center * -2.0, radius),
    ]
    pharmacophore = LigandBasedPharmacophore(None)
    pharmacophore.pharmacophoric_points = points
    return pharmacophore


def test_pharmacophore_len(pharmacophore_three_points):
    assert len(pharmacophore_three_points) == 3


def test_iterate_pharmacophore(pharmacophore_three_points):
    short_names = ["D", "R", "H"]
    ii = 0
    iterated = False
    for point in pharmacophore_three_points:
        assert point.short_name == short_names[ii]
        ii += 1
        iterated = True

    assert iterated


def test_pharmacophore_get_item(pharmacophore_three_points):
    assert pharmacophore_three_points[0].feature_name == "hb donor"
    assert pharmacophore_three_points[1].feature_name == "aromatic ring"
    assert pharmacophore_three_points[2].feature_name == "hydrophobicity"


def test_pharmacophore_equality(pharmacophore_three_points):
    radius = puw.quantity(1.0, "angstroms")
    donor = PharmacophoricPoint("hb donor",
                                puw.quantity([1.0, 1.0, 1.0], "angstroms"),
                                radius)
    pharmacophore_1 = LigandBasedPharmacophore(None)
    pharmacophore_1.pharmacophoric_points = [donor]
    assert not pharmacophore_1 == pharmacophore_three_points

    ring = PharmacophoricPoint("aromatic ring",
                               puw.quantity([2.0, 2.0, 2.0], "angstroms"),
                               radius)
    hydrophobic = PharmacophoricPoint("hydrophobicity",
                                      puw.quantity([-1.0, 2.0, 2.0], "angstroms"),
                                      radius)
    pharmacophore_2 = LigandBasedPharmacophore(None)
    pharmacophore_2.pharmacophoric_points = [donor, ring, hydrophobic]
    assert not pharmacophore_2 == pharmacophore_three_points

    hydrophobic.center = puw.quantity([-2.0, -2.0, -2.0], "angstroms")
    pharmacophore_3 = LigandBasedPharmacophore(None)
    pharmacophore_3.pharmacophoric_points = [donor, ring, hydrophobic]
    assert pharmacophore_3 == pharmacophore_three_points


def test_to_rdkit(pharmacophore_three_points):
    rdkit_ph, radii = pharmacophore_three_points.to_rdkit()

    assert len(radii) == 3
    assert np.allclose(np.array(radii), np.array([1.0, 1.0, 1.0]))

    feats = rdkit_ph.getFeatures()
    assert len(feats) == 3

    acceptor = feats[0]
    assert acceptor.GetFamily() == "Donor"
    assert np.allclose(acceptor.GetPos().x, 1.0)
    assert np.allclose(acceptor.GetPos().y, 1.0)
    assert np.allclose(acceptor.GetPos().z, 1.0)

    ring_1 = feats[1]
    assert ring_1.GetFamily() == "Aromatic"
    assert np.allclose(ring_1.GetPos().x, 2.0)
    assert np.allclose(ring_1.GetPos().y, 2.0)
    assert np.allclose(ring_1.GetPos().z, 2.0)

    ring_2 = feats[2]
    assert ring_2.GetFamily() == "Hydrophobe"
    assert np.allclose(ring_2.GetPos().x, -2.0)
    assert np.allclose(ring_2.GetPos().y, -2.0)
    assert np.allclose(ring_2.GetPos().z, -2.0)


def test_distance_matrix():
    # Test pharmacophore with zero elements
    pharmacophore = LigandBasedPharmacophore(None)
    matrix = pharmacophore.distance_matrix()
    assert matrix.shape == (0, 0)

    radius = puw.quantity(1.0, "angstroms")
    points = [
        PharmacophoricPoint(feat_type="hb acceptor",
                            center=puw.quantity([1, 0, -5], "angstroms"),
                            radius=radius),
        PharmacophoricPoint(feat_type="hb donor",
                            center=puw.quantity([2, 1, 0], "angstroms"),
                            radius=radius),
        PharmacophoricPoint(feat_type="aromatic ring",
                            center=puw.quantity([-3, 2, -1], "angstroms"),
                            radius=radius),
    ]

    pharmacophore = LigandBasedPharmacophore(None)
    pharmacophore.pharmacophoric_points = points
    distance_matrix = pharmacophore.distance_matrix()

    assert distance_matrix.shape == (3, 3)
    assert np.allclose(distance_matrix,
                       np.array([[0, np.sqrt(27), 6],
                                 [np.sqrt(27), 0, np.sqrt(27)],
                                 [6, np.sqrt(27), 0]])
                       )

    points = [
        PharmacophoricPoint(feat_type="hb acceptor",
                            center=puw.quantity([1, 2, 3], "angstroms"),
                            radius=radius),
        PharmacophoricPoint(feat_type="negative charge",
                            center=puw.quantity([0, 0, 0], "angstroms"),
                            radius=radius),
        PharmacophoricPoint(feat_type="positive charge",
                            center=puw.quantity([-1, 0, -1], "angstroms"),
                            radius=radius),
        PharmacophoricPoint(feat_type="aromatic ring",
                            center=puw.quantity([2, -1, 1], "angstroms"),
                            radius=radius),
    ]

    sq = np.sqrt
    pharmacophore = LigandBasedPharmacophore(None)
    pharmacophore.pharmacophoric_points = points
    distance_matrix = pharmacophore.distance_matrix()
    assert distance_matrix.shape == (4, 4)
    assert np.allclose(distance_matrix,
                       np.array([[0, sq(14), sq(24), sq(14)],
                                 [sq(14), 0, sq(2), sq(6)],
                                 [sq(24), sq(2), 0, sq(14)],
                                 [sq(14), sq(6), sq(14), 0]]))


@pytest.mark.skip(reason="Not implemented yet")
def test_to_ligand_scout():
    assert False, "Implement me!"


@pytest.mark.skip(reason="Not implemented yet")
def test_to_moe():
    assert False, "Implement me!"


@pytest.mark.skip(reason="Not implemented yet")
def test_to_pharmagist():
    assert False, "Implement me!"


def test_pharmacophore_string_representation(pharmacophore_three_points):
    assert pharmacophore_three_points.__repr__() == "LigandBasedPharmacophore(n_pharmacophoric_points: 3)"
