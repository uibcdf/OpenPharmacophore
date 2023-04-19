import pyunitwizard as puw
from openpharmacophore import LigandBasedPharmacophore, LigandReceptorPharmacophore, \
    Pharmacophore, PharmacophoricPoint, VirtualScreening


def test_accepts_different_types_of_pharmacophores():

    radius = puw.quantity(1.0, "angstroms")
    points = [
        PharmacophoricPoint("hb donor", puw.quantity([1.0, 1.0, 1.0], "angstroms"), radius),
        PharmacophoricPoint("hb acceptor", puw.quantity([1.0, 1.0, 1.0], "angstroms"), radius),
    ]
    lbp = LigandBasedPharmacophore([])
    lbp._pharmacophores = [Pharmacophore(points)]
    lrp = LigandReceptorPharmacophore(None, None)
    lrp._pharmacophores = [Pharmacophore(points)]
    pharmacophore = Pharmacophore(points)
    screener = VirtualScreening([lbp, lrp, pharmacophore])
    assert screener.n_pharmacophores == 3
    assert screener.pharmacophores[0] == lbp[0]
    assert screener.pharmacophores[1] == lrp[0]
    assert screener.pharmacophores[2] == pharmacophore
