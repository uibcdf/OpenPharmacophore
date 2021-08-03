from openpharmacophore.ligand_based import LigandBasedPharmacophore
from openpharmacophore import pharmacophoric_elements
import pyunitwizard as puw

def test_add_element():
    hb_acceptor = pharmacophoric_elements.HBAcceptorSphere(
        center=puw.quantity([1,0,0], "angstroms"),
        radius=puw.quantity(1.0, "angstroms")
    )
    ring = pharmacophoric_elements.AromaticRingSphere(
        center=puw.quantity([2, 1, 4], "angstroms"),
        radius=puw.quantity(1.0, "angstroms")
    )
    pharmacophore = LigandBasedPharmacophore(elements=[hb_acceptor, ring])
    hydrophobic = pharmacophoric_elements.HydrophobicSphere(
        center=puw.quantity([-1, 0, 2], "angstroms"),
        radius=puw.quantity(1.0, "angstroms")
    )
    pharmacophore.add_element(hydrophobic)
    assert pharmacophore.n_elements == 3
    assert isinstance(pharmacophore.elements[2], pharmacophoric_elements.HydrophobicSphere)

def test_reset():

    hb_acceptor = pharmacophoric_elements.HBAcceptorSphere(
        center=puw.quantity([1,0,0], "angstroms"),
        radius=puw.quantity(1.0, "angstroms")
    )
    ring = pharmacophoric_elements.AromaticRingSphere(
        center=puw.quantity([2, 1, 4], "angstroms"),
        radius=puw.quantity(1.0, "angstroms")
    )
    pharmacophore = LigandBasedPharmacophore(elements=[hb_acceptor, ring])
    pharmacophore._reset()
    assert pharmacophore.n_elements == 0
    assert len(pharmacophore.elements) == 0
    assert pharmacophore.extractor is None
    assert pharmacophore.molecular_system is None

def test_add_to_NGLView():
    pass

def test_show():
    pass