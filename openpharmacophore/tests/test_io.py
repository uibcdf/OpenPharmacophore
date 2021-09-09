from openpharmacophore.io.mol2 import load_mol2_file
from openpharmacophore.io.moe import from_moe, to_moe
from openpharmacophore.io.ligandscout import from_ligandscout, to_ligandscout
from openpharmacophore.io.pharmagist import read_pharmagist, to_pharmagist
from openpharmacophore.io.pharmer import from_pharmer, to_pharmer
import openpharmacophore.pharmacophoric_elements as phe
from openpharmacophore.ligand_based import LigandBasedPharmacophore
from openpharmacophore.structured_based import StructuredBasedPharmacophore

import numpy as np
import pyunitwizard as puw
import pytest

import datetime
import os

@pytest.fixture
def two_element_pharmacophore():
    """Returns a pharmacophore with an aromatic ring and an hb acceptor"""
    radius = puw.quantity(1.0, "angstroms")
    ring = phe.AromaticRingSphereAndVector(puw.quantity([1,0,0], "angstroms"), 
                                            radius, [0, 0, 1])
    acceptor = phe.HBAcceptorSphere(puw.quantity([1,2,2], "angstroms"), 
                                    radius)
    pharmacophore = StructuredBasedPharmacophore(elements=[ring, acceptor])
    return pharmacophore

@pytest.fixture
def three_element_pharmacophore():
    radius = puw.quantity(1.0, "angstroms")
    ring = phe.AromaticRingSphereAndVector(puw.quantity([1,0,0], "angstroms"), 
                                            radius, [0, 0, 1])
    acceptor = phe.HBAcceptorSphereAndVector(puw.quantity([1,2,2], "angstroms"), 
                                    radius, [0,1,1])
    excluded = phe.ExcludedVolumeSphere(puw.quantity([2,1,2], "angstroms"), radius)
    pharmacophore = StructuredBasedPharmacophore(elements=[ring, acceptor, excluded])
    return pharmacophore

def test_from_pharmer():

    points, molecular_system = from_pharmer("./openpharmacophore/data/pharmacophores/pharmer/pharmer.json", 
                                            load_mol_sys=False)
    assert len(points) == 19
    assert molecular_system is None
    assert isinstance(points[0], 
                        phe.aromatic_ring.AromaticRingSphereAndVector)

def test_to_pharmer(two_element_pharmacophore):
    
    pharmer = to_pharmer(two_element_pharmacophore, "temp.json", testing=True)

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

def test_load_mol2_file():
    fname = "./openpharmacophore/data/ligands/ace.mol2"
    molecules = load_mol2_file(fname=fname)

    assert len(molecules) == 3
    assert molecules[0].GetNumAtoms() == 14
    assert molecules[1].GetNumAtoms() == 25
    assert molecules[2].GetNumAtoms() == 29

@pytest.mark.parametrize("fname,index",[
    ("elastase.mol2", None),
    ("elastase.mol2", 0),
    ("streptadivin.mol2", None)
])
def test_read_pharmagist(fname, index):
    path = "./openpharmacophore/data/pharmacophores/pharmagist"
    fpath = os.path.join(path, fname)
    
    result = read_pharmagist(fpath, index)
    if index is None:
        assert isinstance(result[0], LigandBasedPharmacophore)
        if fname == "elastase.mol2":
            assert result[0].n_elements == 4
            assert len(result) == 8
        else:
            assert result[0].n_elements == 9
            assert len(result) == 6
    else:
        assert len(result) == 4
        assert isinstance(result[0], phe.HBAcceptorSphere)

def test_to_pharmagist(three_element_pharmacophore):
    mol2_list = to_pharmagist(three_element_pharmacophore, file_name=None, testing=True)
    expected_output = ['@<TRIPOS>MOLECULE\n',
                        '@<TRIPOS>ATOM\n',
                        '      1 AR           1.0000    0.0000    0.0000   AR     0   AR      0.0000\n',
                        '      2 ACC          1.0000    2.0000    2.0000   HB     1   HB      0.0000\n',
                        '@<TRIPOS>BOND\n']
    assert mol2_list == expected_output

def test_from_ligandscout():
    points = from_ligandscout("./openpharmacophore/data/pharmacophores/ligandscout/pharmacophore.pml")
    assert len(points) == 4

    neg_ion = points[0]  
    assert isinstance(neg_ion, phe.NegativeChargeSphere)
    assert np.all(
            np.around(puw.get_value(neg_ion.center, "angstroms"), 1) == np.array([-8.0, 10.0, -9.5])
            )
    assert puw.get_value(neg_ion.radius, "angstroms") == 1.5
    
    donor = points[1]
    assert isinstance(donor, phe.HBDonorSphereAndVector)
    assert np.all(
            np.around(puw.get_value(donor.center, "angstroms"), 1) ==  np.array([-8.0, 2.0, -10.0])
            )
    assert puw.get_value(donor.radius, "angstroms") == 1.5
    dir_expected = np.array([-5.690445899963379, 0.5822541117668152, -10.5515718460083])
    dir_expected /= np.linalg.norm(dir_expected)
    assert np.all(
            np.around(donor.direction, 2) == np.around(dir_expected, 2)
            )
    
    ring = points[2]
    assert isinstance(ring, phe.AromaticRingSphereAndVector)
    assert np.all(
            np.around(puw.get_value(ring.center, "angstroms"), 1) ==  np.array([0.0, 6.5, -3.0])
            )
    assert puw.get_value(ring.radius, "angstroms") == 1.5
    dir_expected = np.array([-3.8126893043518066, 1.7578959465026855, 0.6093783378601074])
    dir_expected /= np.linalg.norm(dir_expected)
    assert np.all(
            np.around(ring.direction, 2) == np.around(dir_expected, 2)
            )
    
    excluded_vol = points[3]
    assert isinstance(excluded_vol, phe.ExcludedVolumeSphere)
    assert np.all(
            np.around(puw.get_value(excluded_vol.center, "angstroms"), 1) ==  np.array([5.5, 4.5, -2.0])
            )
    assert round(puw.get_value(excluded_vol.radius, "angstroms"), 1) == 1.0

def test_to_ligandscout(three_element_pharmacophore):
    expected_string = b'<pharmacophore name="pharmacophore.pml" pharmacophoreType="LIGAND_SCOUT"><plane disabled="false" featureId="ai_1" name="AR" optional="false" weight="1.0"><position tolerance="0.9999999999999999" x3="0.9999999999999999" y3="0.0" z3="0.0" /><normal tolerance="0.9999999999999999" x3="0.0" y3="0.0" z3="1.0" /></plane><vector disabled="false" featureId="ha_2" hasSyntheticProjectedPoint="false" name="HBA" optional="false" pointsToLigand="false" weight="1.0"><origin tolerance="0.9999999999999999" x3="0.9999999999999999" y3="1.9999999999999998" z3="1.9999999999999998" /><target tolerance="0.9999999999999999" x3="0.9999999999999999" y3="1.2928932188134523" z3="1.2928932188134523" /></vector><volume disabled="false" featureId="ev_3" optional="false" type="exclusion" weight="1.0"><position tolerance="0.9999999999999999" x3="1.9999999999999998" y3="0.9999999999999999" z3="1.9999999999999998" /></volume></pharmacophore>'
    pml_string = to_ligandscout(three_element_pharmacophore, file_name=None, testing=True)

    assert pml_string == expected_string

def test_from_moe():
    file_name = "./openpharmacophore/data/pharmacophores/moe/gmp.ph4"
    points = from_moe(file_name)
    assert len(points) == 10

    donor = points[0]
    assert isinstance(donor, phe.HBDonorSphere)
    assert np.all(
        np.around(puw.get_value(donor.center, "angstroms"), 2) == np.around(np.array([1.71, 1.43075, -1.4255]), 2)
        )
    assert np.around(puw.get_value(donor.radius, "angstroms"), 2) == 0.51
        
    
    hyd = points[1]
    assert isinstance(hyd, phe.HydrophobicSphere)
    assert np.all(
        np.around(puw.get_value(hyd.center, "angstroms"), 2) == np.around(np.array([2.7895, 2.4035, -1.40875]), 2)
        )
    assert np.around(puw.get_value(hyd.radius, "angstroms"), 2) == 0.55
    
    acceptor = points[2]
    assert isinstance(acceptor, phe.HBAcceptorSphere)
    assert np.all(
        np.around(puw.get_value(acceptor.center, "angstroms"), 2) == np.around(np.array([0.312, 3.0175, -2.44825]), 2)
        )
    assert np.around(puw.get_value(acceptor.radius, "angstroms"), 2) == 0.57
    
    aromatic_1 = points[3]
    assert isinstance(aromatic_1, phe.AromaticRingSphere)
    assert np.all(
        np.around(puw.get_value(aromatic_1.center, "angstroms"), 2) == 
        np.around(np.array([-0.748458333333333, 2.13108333333333, -2.490375]), 2)
        )
    assert np.around(puw.get_value(aromatic_1.radius, "angstroms"), 2) == 0.58
        
    
    aromatic_2 = points[4]
    assert isinstance(aromatic_2, phe.AromaticRingSphere)
    assert np.all(
        np.around(puw.get_value(aromatic_2.center, "angstroms"), 2) == 
        np.around(np.array([-1.719625, -0.0273333333333334, -2.055625]), 2)
        )
    assert np.around(puw.get_value(aromatic_2.radius, "angstroms"), 2) == 0.6

    
    aromatic_3 = points[5]
    assert isinstance(aromatic_3, phe.AromaticRingSphere)
    assert np.all(
        np.around(puw.get_value(aromatic_3.center, "angstroms"), 2) == 
        np.around(np.array([5.20029166666667, 1.25479166666667, -0.199041666666667]), 2)
        )
    assert np.around(puw.get_value(aromatic_3.radius, "angstroms"), 2) == 0.61
    
    acceptor = points[6]
    assert isinstance(acceptor, phe.HBAcceptorSphere)
    assert np.all(
        np.around(puw.get_value(acceptor.center, "angstroms"), 2) == 
        np.around(np.array([-1.95875, 2.536, -3.03625]), 2)
        )
    assert np.around(puw.get_value(acceptor.radius, "angstroms"), 2) == 0.62
        
    
    hyd_2 = points[7]
    assert isinstance(hyd_2, phe.HydrophobicSphere)
    assert np.all(
        np.around(puw.get_value(hyd_2.center, "angstroms"), 2) == np.around(np.array([-1.54725, -2.979375, -0.961875]), 2)
        )
    assert np.around(puw.get_value(hyd_2.radius, "angstroms"), 2) == 0.74
      
    
    acceptor_2 = points[8]
    assert isinstance(acceptor_2, phe.HBAcceptorSphere)
    assert np.all(
        np.around(puw.get_value(acceptor_2.center, "angstroms"), 2) == 
        np.around(np.array([-0.755095833333333, 6.3286375, -3.96758333333333]), 2)
        )
    assert np.around(puw.get_value(acceptor_2.radius, "angstroms"), 2) == 1.25
    
def test_to_moe(three_element_pharmacophore):

    pharmacophore_str = to_moe(three_element_pharmacophore, file_name=None, testing=True)

    now = datetime.datetime.now()
    month = str(now.month)
    year = str(now.year)
    expected_str = f'#moe:ph4que {year}.{month}\n#pharmacophore 5 tag t value *\nscheme t Unified matchsize i 0 title t s $\n#feature 3 expr tt color ix x r y r z r r r ebits ix gbits ix\nAro df2f2 0.9999999999999999 0.0 0.0 0.9999999999999999 0 300 Acc df2f2 0.9999999999999999 1.9999999999999998 1.9999999999999998 0.9999999999999999 0 300 \n#volumesphere 90 x r y r z r r r\n1.9999999999999998 0.9999999999999999 1.9999999999999998 0.9999999999999999 \n#endpharmacophore'

    assert pharmacophore_str == expected_str

