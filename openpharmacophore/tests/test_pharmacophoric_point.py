from openpharmacophore.pharmacophoric_point import PharmacophoricPoint
import numpy as np
import pyunitwizard as puw

def test_PharmacophoricPoint():

    feat_name = "hb donor"
    center = puw.quantity([1.0, 1.0, 1.0], "angstroms")
    radius = puw.quantity(1.0, "angstroms")
    atom_inxs = (3, 4, 5, 6)

    donor_1 = PharmacophoricPoint(feat_name, center, radius, None, atom_inxs)
    donor_2 = PharmacophoricPoint(feat_name, center, radius, 
                direction=np.array([1.0, 1.0, 1.0]), atoms_inxs=atom_inxs)

    assert donor_1.has_direction == False
    assert donor_1.element_name == "HbDonorSphere"
    assert np.allclose(donor_1.get_center(), np.array([1.0, 1.0, 1.0]))
    assert np.allclose(donor_1.get_radius(), 1.0)
    assert atom_inxs == (3, 4, 5, 6) 

    feat_name = "aromatic ring"
    center = puw.quantity([1.5, -2.0, 3.2], "angstroms")
    radius = puw.quantity(1.5, "angstroms")
    direction = np.array([1.0, 1.0, 1.0])
    atom_inxs = None

    ring = PharmacophoricPoint(feat_name, center, radius, direction, atom_inxs)
    assert ring.has_direction == True
    assert ring.element_name == "AromaticRingSphereAndVector"
    assert ring.atoms_inxs is None
    assert np.allclose(ring.get_center(), np.array([1.5, -2.0, 3.2]))
    assert np.allclose(ring.get_radius(), 1.5) 
    assert np.allclose(ring.get_direction(), 
            np.array([[1.0, 1.0, 1.0]]) / np.linalg.norm(np.array([1.0, 1.0, 1.0])))

    assert donor_1 != ring
    assert donor_2 != donor_1