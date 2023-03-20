import numpy as np
import pyunitwizard as puw

from openpharmacophore.molecular_systems.chem_feats import ChemFeatContainer, ChemFeat, HBDonor, AromaticRing
from openpharmacophore import LigandReceptorPharmacophore as LRP
from openpharmacophore import PharmacophoricPoint


def test_hbond_angle():
    donor_acc_dist = puw.quantity(0.5, "angstroms")
    don_hyd_dist = puw.quantity(0.3535533906, "angstroms")
    acc_hyd_dist = puw.quantity(0.3535533906, "angstroms")

    angle = LRP._hbond_angle(donor_acc_dist, don_hyd_dist, acc_hyd_dist)
    assert np.allclose(puw.get_value(angle), 90)


def test_hbond_pharmacophoric_points_distance_cutoff_exceeded():
    donors = [HBDonor(coords=puw.quantity(np.array([1., 0., 3.]), "angstroms"),
                      hyd=puw.quantity(np.array([0., 0., 0.]), "angstroms"))]
    acceptors = [
        ChemFeat(type="hb acceptor", coords=puw.quantity(np.array([10., 10., 10.]), "angstroms")),
    ]
    points = LRP._hbond_pharmacophoric_points(donors, acceptors, True)
    assert len(points) == 0


def test_hbond_pharmacophoric_points_angle_below_cutoff():
    lig_feats = [HBDonor(coords=puw.quantity(np.array([0., 0., 0.]), "angstroms"),
                         hyd=puw.quantity(np.array([0.25, 0.25, 0.]), "angstroms"))]
    rec_feats = [
        ChemFeat(type="hb acceptor", coords=puw.quantity(np.array([0.3, 0., 0.]), "angstroms")),
    ]
    # D-H--A angle should be about 56°
    points = LRP._hbond_pharmacophoric_points(lig_feats, rec_feats, True)
    assert len(points) == 0


def test_hbond_pharmacophoric_points():
    lig_feats = [HBDonor(coords=puw.quantity(np.array([0., 0., 0.]), "angstroms"),
                         hyd=puw.quantity(np.array([0.25, 0.5, 0.]), "angstroms"))]
    rec_feats = [
        ChemFeat(type="hb acceptor", coords=puw.quantity(np.array([0.75, 0.75, 0.]), "angstroms")),
    ]

    points = LRP._hbond_pharmacophoric_points(lig_feats, rec_feats, donors_in_ligand=True)
    assert len(points) == 1
    assert points[0].feature_name == "hb donor"
    assert np.all(puw.get_value(points[0].center) == np.array([0., 0., 0.]))
    assert np.allclose(points[0].direction, np.array([0.89443, 0.44721, 0.]))


def test_aromatic_pharmacophoric_points_exceeds_max_distance():
    ligand_feats = [AromaticRing(puw.quantity(np.array([0., 0., 0.]), "angstroms"))]
    receptor_feats = [AromaticRing(puw.quantity(np.array([10., 10., 10.]), "angstroms"))]
    points = LRP._aromatic_pharmacophoric_points(ligand_feats, receptor_feats)
    assert len(points) == 0


def aromatic_ligand_feats():
    return [
        AromaticRing(
            coords=puw.quantity(np.array([0., 0., 0.]), "angstroms"),
            normal=np.array([0., 0., 1.]))
    ]


def aromatic_receptor_feats():
    return [
        AromaticRing(
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
        AromaticRing(
            coords=puw.quantity(np.array([5., 0., 3.]), "angstroms"),
            normal=np.array([0., 0., 1.]))
    ]

    points = LRP._aromatic_pharmacophoric_points(lig_feats, rec_feats)
    assert len(points) == 0


def test_aromatic_pharmacophoric_points_t_stack():
    lig_feats = aromatic_ligand_feats()
    rec_feats = [
        AromaticRing(
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
        ChemFeat(type="positive charge", coords=puw.quantity(np.array([0.] * 3), "angstroms")),
        ChemFeat(type="positive charge", coords=puw.quantity(np.array([-7.] * 3), "angstroms"))
    ]
    receptor_feats = [
        ChemFeat(type="negative charge", coords=puw.quantity(np.array([7.] * 3), "angstroms")),
        ChemFeat(type="negative charge", coords=puw.quantity(np.array([2.] * 3), "angstroms"))
    ]

    points = LRP._points_from_distance_rule(
        ligand_feats, receptor_feats, "positive charge", LRP.CHARGE_DIST_MAX
    )

    assert len(points) == 1
    assert np.all(puw.get_value(points[0].center) == np.array([0.] * 3))


class FakeLigand:

    def get_chem_feats_with_directionality(self, conf_ind):
        donor = HBDonor(coords=puw.quantity(np.array([0., 0., 0.]), "angstroms"),
                        hyd=puw.quantity(np.array([0.25, 0.5, 0.]), "angstroms"))
        ring = AromaticRing(coords=puw.quantity(np.array([0., 0., 0.]), "angstroms"),
                            normal=np.array([0., 0., 1.]))
        pos = ChemFeat(type="positive charge", coords=puw.quantity(np.array([0.] * 3), "angstroms"))
        neg = ChemFeat(type="negative charge", coords=puw.quantity(np.array([0.] * 3), "angstroms"))
        hyd = ChemFeat(type="hydrophobicity", coords=puw.quantity(np.array([0.] * 3), "angstroms"))

        container = ChemFeatContainer()
        container.add_feats([donor, ring, pos, neg, hyd])
        return container


class FakeBindingSite:

    def get_chem_feats(self, frame):
        acceptor = ChemFeat(type="hb acceptor", coords=puw.quantity(np.array([0.75, 0.75, 0.]), "angstroms"))
        ring = AromaticRing(
            coords=puw.quantity(np.array([1., 0., 3.]), "angstroms"),
            normal=np.array([0., 0., 1.]))
        pos = ChemFeat(type="positive charge", coords=puw.quantity(np.array([2.] * 3), "angstroms"))
        neg = ChemFeat(type="negative charge", coords=puw.quantity(np.array([2.] * 3), "angstroms"))
        hyd = ChemFeat(type="hydrophobicity", coords=puw.quantity(np.array([2.] * 3), "angstroms"))

        container = ChemFeatContainer()
        container.add_feats([acceptor, ring, pos, neg, hyd])
        return container


def test_extract():
    pharmacophore = LRP(FakeBindingSite(), FakeLigand())
    pharmacophore.extract()

    assert len(pharmacophore) == 1
    assert len(pharmacophore[0]) == 5

    assert pharmacophore[0].ref_struct == 0
    assert pharmacophore[0].ref_mol == 0
