import openpharmacophore as oph
import pyunitwizard as puw
from openpharmacophore.io.pharmer import from_pharmer, to_pharmer
from openpharmacophore.structured_based import StructuredBasedPharmacophore

def test_from_pharmer():

    pharmacophore = from_pharmer("./openpharmacophore/data/pharmer.json", load_mol_sys=False)
    assert pharmacophore.n_elements == 19
    assert pharmacophore.molecular_system is None
    assert isinstance(pharmacophore.elements[0], 
                        oph.pharmacophoric_elements.aromatic_ring.AromaticRingSphereAndVector)

def test_to_pharmer():
    
    # Create pharmacohore with two points
    pharmacophore = StructuredBasedPharmacophore()
    radius = puw.quantity(1.0, "angstroms")
    ring_center = puw.quantity([1,0,0], "angstroms")
    ring_direction = [0,0,1]
    hb_center = puw.quantity([1,2,2], "angstroms")
    pharmacophore.add_element(oph.pharmacophoric_elements.AromaticRingSphereAndVector(ring_center, radius, ring_direction))
    pharmacophore.add_element(oph.pharmacophoric_elements.HBAcceptorSphere(hb_center, radius))

    pharmer = to_pharmer(pharmacophore, "temp.json", return_dict=True)

    # Expected output from to_pharmer
    expected = {}
    expected["points"] = []

    aromatic = {}
    aromatic["name"] = "Aromatic"
    aromatic["hasvec"] = True
    aromatic["svector"] = {}
    aromatic["svector"]["x"] = 0.0
    aromatic["svector"]["y"] = 0.0
    aromatic["svector"]["z"] = 1.0
    aromatic["x"] = 0.9999999999999999
    aromatic["y"] = 0.0
    aromatic["z"] = 0.0
    aromatic["radius"] = 0.9999999999999999
    aromatic["enabled"] = True
    aromatic["vector_on"] = 0
    aromatic["minsize"] = ""
    aromatic["maxsize"] = ""
    aromatic["selected"] = False

    acceptor = {}
    acceptor["name"] = "HydrogenAcceptor"
    acceptor["hasvec"] = False
    acceptor["svector"] = {}
    acceptor["svector"]["x"] = 1
    acceptor["svector"]["y"] = 0
    acceptor["svector"]["z"] = 0
    acceptor["x"] = 0.9999999999999999
    acceptor["y"] = 1.9999999999999998
    acceptor["z"] = 1.9999999999999998
    acceptor["radius"] = 0.9999999999999999
    acceptor["enabled"] = True
    acceptor["vector_on"] = 0
    acceptor["minsize"] = ""
    acceptor["maxsize"] = ""
    acceptor["selected"] = False

    expected["points"].append(aromatic) 
    expected["points"].append(acceptor)

    assert pharmer == expected