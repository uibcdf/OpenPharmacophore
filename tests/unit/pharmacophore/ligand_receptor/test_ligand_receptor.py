import numpy as np
import pyunitwizard as puw
import pytest

from openpharmacophore.molecular_systems.chem_feats import ChemFeat
from openpharmacophore import LigandReceptorPharmacophore as LRP
from openpharmacophore import PharmacophoricPoint


@pytest.mark.skip(reason="temp")
def test_hb_donor_pharmacophoric_points():
    pass


@pytest.mark.skip(reason="temp")
def test_hb_acceptor_pharmacophoric_points():
    pass


def test_aromatic_pharmacophoric_points_exceeds_max_distance():
    ligand_feats = [ChemFeat("aromatic ring", puw.quantity(np.array([0., 0., 0.]), "angstroms"))]
    receptor_feats = [ChemFeat("aromatic ring", puw.quantity(np.array([10., 10., 10.]), "angstroms"))]
    points = LRP._aromatic_pharmacophoric_points(ligand_feats, receptor_feats)
    assert len(points) == 0


def aromatic_ligand_feats():
    return [
        ChemFeat(
            type="aromatic ring",
            coords=puw.quantity(np.array([0., 0., 0.]), "angstroms"),
            normal=np.array([0., 0., 1.]))
    ]


def aromatic_receptor_feats():
    return [
        ChemFeat(
            type="aromatic ring",
            coords=puw.quantity(np.array([1., 0., 3.]), "angstroms"),
            normal=np.array([0., 0., 1.]))
    ]


def test_aromatic_pharmacophoric_points_p_stack():
    points = LRP._aromatic_pharmacophoric_points(
        aromatic_ligand_feats(), aromatic_receptor_feats()
    )
    assert len(points) == 1
    assert points[0].feature_name == "aromatic ring"

    assert np.all(puw.get_value(points[0].center) == np.array([0.] * 3))
    assert points[0].has_direction

    assert np.allclose(points[0].direction,
                       np.array([1 / np.sqrt(10), 0, 3 / np.sqrt(10)]))


def test_aromatic_pharmacophoric_points_max_offset_exceeded():
    lig_feats = aromatic_ligand_feats()
    rec_feats = [
        ChemFeat(type="aromatic ring",
                 coords=puw.quantity(np.array([5., 0., 3.]), "angstroms"),
                 normal=np.array([0., 0., 1.]))
    ]

    points = LRP._aromatic_pharmacophoric_points(lig_feats, rec_feats)
    assert len(points) == 0


def test_aromatic_pharmacophoric_points_t_stack():
    lig_feats = aromatic_ligand_feats()
    rec_feats = [
        ChemFeat(type="aromatic ring",
                 coords=puw.quantity(np.array([2., 0., 0.]), "angstroms"),
                 normal=np.array([1., 0., 0.]))
    ]

    points = LRP._aromatic_pharmacophoric_points(lig_feats, rec_feats)
    assert len(points) == 1
    assert points[0].feature_name == "aromatic ring"

    assert np.all(puw.get_value(points[0].center) == np.array([0.] * 3))
    assert points[0].has_direction

    assert np.allclose(points[0].direction, np.array([1., 0., 0.]))


def test_merge_hydrophobic_points():
    radius = puw.quantity(1.0, "angstroms")
    centers = [
        puw.quantity(np.array([0., 0., 0.]), "angstroms"),
        puw.quantity(np.array([1., 1., 0.]), "angstroms"),
        puw.quantity(np.array([0., 1., 0.]), "angstroms"),
        puw.quantity(np.array([0., 2.5, 0.]), "angstroms"),
    ]
    points = [
        PharmacophoricPoint("hydrophobicity", c, radius) for c in centers
    ]
    merged = LRP._merge_hydrophobic_points(points, radius)
    assert len(merged) == 2
    assert all([p.feature_name == "hydrophobicity" for p in merged])
    assert np.allclose(puw.get_value(merged[0].center), np.array([1 / 3, 2 / 3, 0]))
    assert np.allclose(puw.get_value(merged[1].center), np.array([1 / 3, 3 / 2, 0]))


def test_pharmacophoric_points_from_distance_rule():
    ligand_feats = [
        ChemFeat("positive charge", puw.quantity(np.array([0.] * 3), "angstroms")),
        ChemFeat("positive charge", puw.quantity(np.array([-7.] * 3), "angstroms"))
    ]
    receptor_feats = [
        ChemFeat("negative charge", puw.quantity(np.array([7.] * 3), "angstroms")),
        ChemFeat("negative charge", puw.quantity(np.array([2.] * 3), "angstroms"))
    ]

    points = LRP._points_from_distance_rule(
        ligand_feats, receptor_feats, "positive charge", LRP.CHARGE_DIST_MAX
    )

    assert len(points) == 1
    assert np.all(puw.get_value(points[0].center) == np.array([0.] * 3))
