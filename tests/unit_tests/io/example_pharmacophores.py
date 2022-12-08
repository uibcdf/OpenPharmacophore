from openpharmacophore import PharmacophoricPoint
import pyunitwizard as puw


def two_element_pharmacophore():
    """Returns a pharmacophore with an aromatic ring and a hb acceptor"""
    radius = puw.quantity(1.0, "angstroms")
    ring = PharmacophoricPoint(
        feat_type="aromatic ring",
        center=puw.quantity([1, 0, 0], "angstroms"),
        radius=radius,
        direction=[0, 0, 1],
    )
    acceptor = PharmacophoricPoint(
        feat_type="hb acceptor",
        center=puw.quantity([1, 2, 2], "angstroms"),
        direction=[0, 1, 1],
        radius=radius)
    return [ring, acceptor]


def three_element_pharmacophore():
    radius = puw.quantity(1.0, "angstroms")
    ring = PharmacophoricPoint(
        feat_type="aromatic ring",
        center=puw.quantity([1, 0, 0], "angstroms"),
        radius=radius,
        direction=[0, 0, 1],
    )
    acceptor = PharmacophoricPoint(
        feat_type="hb acceptor",
        center=puw.quantity([1, 2, 2], "angstroms"),
        radius=radius,
        direction=[0, 1, 1]
    )
    excluded = PharmacophoricPoint(
        feat_type="excluded volume",
        center=puw.quantity([2, 1, 2], "angstroms"),
        radius=radius)
    return [acceptor, excluded, ring]


def five_element_pharmacophore():
    radius = puw.quantity(1.0, "angstroms")
    ring_1 = PharmacophoricPoint(
        feat_type="aromatic ring",
        center=puw.quantity([1, 0, 0], "angstroms"),
        radius=radius,
        direction=[0, 0, 1],
    )
    ring_2 = PharmacophoricPoint(
        feat_type="aromatic ring",
        center=puw.quantity([0, 1, 2], "angstroms"),
        radius=puw.quantity(1.0, "angstroms"),
        direction=[0, 0, 1],
    )
    acceptor = PharmacophoricPoint(
        feat_type="hb acceptor",
        center=puw.quantity([1, 2, 2], "angstroms"),
        radius=radius,
        direction=[0, 1, 1],
    )
    hb_donor = PharmacophoricPoint(
        feat_type="hb donor",
        center=puw.quantity([1, 2, 2], "angstroms"),
        radius=radius,
        direction=[0, 1, 1]
    )
    hydrophobicity = PharmacophoricPoint(
        feat_type="hydrophobicity",
        center=puw.quantity([-1, 2, 2], "angstroms"),
        radius=radius,
    )
    return [acceptor, hb_donor,
            hydrophobicity, ring_1, ring_2]
