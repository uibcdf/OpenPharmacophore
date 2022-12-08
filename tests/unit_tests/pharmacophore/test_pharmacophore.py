from openpharmacophore.pharmacophore.pharmacophore import Pharmacophore
from openpharmacophore import PharmacophoricPoint
import pyunitwizard as puw


def test_init_pharmacphore():
    radius = puw.quantity(1.0, "angstroms")
    donor = PharmacophoricPoint(
        "hb donor",
        puw.quantity([1.0, 1.0, 1.0], "angstroms"),
        radius
    )
    pharma = Pharmacophore([donor], score=0.5, ref_mol=0, ref_struct=2)
    assert pharma.score == 0.5
    assert pharma.ref_mol == 0
    assert pharma.ref_struct == 2


def test_add_points_to_pharmacophore():
    radius = puw.quantity(1.0, "angstroms")
    donor = PharmacophoricPoint(
        "hb donor",
        puw.quantity([1.0, 1.0, 1.0], "angstroms"),
        radius
    )
    ring = PharmacophoricPoint(
        "aromatic ring",
        puw.quantity([1.0, 1.0, 1.0], "angstroms"),
        radius
    )
    pharma = Pharmacophore([donor])
    pharma.add(ring)

    assert len(pharma) == 2
    assert pharma[0].feature_name == "hb donor"
    assert pharma[1].feature_name == "aromatic ring"


def test_remove_point_from_pharmacophore():
    radius = puw.quantity(1.0, "angstroms")
    donor = PharmacophoricPoint(
        "hb donor",
        puw.quantity([1.0, 1.0, 1.0], "angstroms"),
        radius
    )
    ring = PharmacophoricPoint(
        "aromatic ring",
        puw.quantity([1.0, 1.0, 1.0], "angstroms"),
        radius
    )
    pharma = Pharmacophore([donor, ring])

    pharma.remove(0)
    assert len(pharma) == 1
    assert pharma[0].feature_name == "aromatic ring"
