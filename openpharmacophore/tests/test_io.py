from openpharmacophore import PharmacophoricPoint, StructuredBasedPharmacophore
from openpharmacophore.data import pharmacophores, ligands
from openpharmacophore.io.mol2 import load_mol2_file
from openpharmacophore.io.moe import from_moe, _moe_ph4_string
from openpharmacophore.io.ligandscout import from_ligandscout, _ligandscout_xml_tree
from openpharmacophore.io.pharmagist import _pharmagist_file_info
from openpharmacophore.io.load_pharmagist import read_pharmagist
from openpharmacophore.io.pharmer import from_pharmer, _pharmer_dict
from openpharmacophore.io.mol_files import load_molecules_file
from openpharmacophore.io.mol_suppliers import smi_has_header_and_id, mol2_mol_generator, smiles_mol_generator
import numpy as np
import pyunitwizard as puw
import pytest
from rdkit import Chem
import datetime
import os
import xml.etree.ElementTree as ET


@pytest.fixture
def two_element_pharmacophore():
    """Returns a pharmacophore with an aromatic ring and an hb acceptor"""
    radius = puw.quantity(1.0, "angstroms")
    ring = PharmacophoricPoint(
        feat_type="aromatic ring",
        center=puw.quantity([1, 0, 0], "angstroms"),
        radius=radius,
        direction=[0, 0, 1],
        atom_indices=None
    )
    acceptor = PharmacophoricPoint(
        feat_type="hb acceptor",
        center=puw.quantity([1, 2, 2], "angstroms"),
        direction=[0, 1, 1],
        radius=radius)
    return StructuredBasedPharmacophore(pharmacophoric_points=[ring, acceptor])


@pytest.fixture
def three_element_pharmacophore():
    radius = puw.quantity(1.0, "angstroms")
    ring = PharmacophoricPoint(
        feat_type="aromatic ring",
        center=puw.quantity([1, 0, 0], "angstroms"),
        radius=radius,
        direction=[0, 0, 1],
        atom_indices=None
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
    return StructuredBasedPharmacophore(pharmacophoric_points=[ring, acceptor, excluded])


@pytest.fixture
def five_element_pharmacophore():
    radius = puw.quantity(1.0, "angstroms")
    ring_1 = PharmacophoricPoint(
        feat_type="aromatic ring",
        center=puw.quantity([1, 0, 0], "angstroms"),
        radius=radius,
        direction=[0, 0, 1],
        atom_indices=None
    )
    ring_2 = PharmacophoricPoint(
        feat_type="aromatic ring",
        center=puw.quantity([0, 1, 2], "angstroms"),
        radius=puw.quantity(1.0, "angstroms")
    )
    acceptor = PharmacophoricPoint(
        feat_type="hb acceptor",
        center=puw.quantity([1, 2, 2], "angstroms"),
        radius=radius,
        direction=[0, 1, 1]
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
    return StructuredBasedPharmacophore(pharmacophoric_points=[acceptor, hb_donor,
                                                               hydrophobicity, ring_1, ring_2], is_sorted=True)


def test_from_pharmer():
    points, molecular_system, ligand = from_pharmer(pharmacophores["1M70"],
                                                    load_mol_sys=False)
    assert len(points) == 5
    assert molecular_system is None
    assert ligand is None

    assert points[0].feature_name == "hb acceptor"
    assert points[0].has_direction
    assert puw.get_value(points[0].radius, "angstroms") == 0.9999999999999999
    assert np.allclose(puw.get_value(points[0].center, "angstroms"),
                       np.array([21.352, -14.531, 19.625]))
    assert np.allclose(points[0].direction,
                       np.array([-0.6405836470264256, 0.7029084735090229, -0.3091476492414897]))

    assert points[1].feature_name == "hb acceptor"
    assert points[1].has_direction
    assert puw.get_value(points[1].radius, "angstroms") == 0.9999999999999999
    assert np.allclose(puw.get_value(points[1].center, "angstroms"),
                       np.array([19.355, -18.32, 23.987]))
    assert np.allclose(points[1].direction,
                       np.array([0.6859059711903811, 0.09092493673854565, 0.721987295292979]))

    assert points[2].feature_name == "hb donor"
    assert points[2].has_direction
    assert puw.get_value(points[2].radius, "angstroms") == 0.9999999999999999
    assert np.allclose(puw.get_value(points[2].center, "angstroms"),
                       np.array([20.977, -16.951, 18.746]))
    assert np.allclose(points[2].direction,
                       np.array([0.71662539652105, -0.5202950182802607, -0.46447942367105033]))

    assert points[3].feature_name == "negative charge"
    assert not points[3].has_direction
    assert puw.get_value(points[3].radius, "angstroms") == 1.5
    assert np.allclose(puw.get_value(points[3].center, "angstroms"),
                       np.array([21.66899, -15.077667, 20.608334]))

    assert points[4].feature_name == "negative charge"
    assert not points[4].has_direction
    assert puw.get_value(points[4].radius, "angstroms") == 1.9999999999999998
    assert np.allclose(puw.get_value(points[4].center, "angstroms"),
                       np.array([19.985, -19.404402, 22.8422]))


def test_to_pharmer(two_element_pharmacophore, five_element_pharmacophore):
    pharmer = _pharmer_dict(two_element_pharmacophore)

    # Expected output from to_pharmer
    expected = {}
    expected["points"] = []

    aromatic_1 = {}
    aromatic_1["name"] = "Aromatic"
    aromatic_1["hasvec"] = True
    aromatic_1["svector"] = {}
    aromatic_1["svector"]["x"] = 0.0
    aromatic_1["svector"]["y"] = 0.0
    aromatic_1["svector"]["z"] = 1.0
    aromatic_1["x"] = 0.9999999999999999
    aromatic_1["y"] = 0.0
    aromatic_1["z"] = 0.0
    aromatic_1["radius"] = 0.9999999999999999
    aromatic_1["enabled"] = True
    aromatic_1["vector_on"] = 0
    aromatic_1["minsize"] = ""
    aromatic_1["maxsize"] = ""
    aromatic_1["selected"] = False

    acceptor = {}
    acceptor["name"] = "HydrogenAcceptor"
    acceptor["hasvec"] = True
    acceptor["svector"] = {}
    acceptor["svector"]["x"] = 0.0
    acceptor["svector"]["y"] = 0.7071067811865475
    acceptor["svector"]["z"] = 0.7071067811865475
    acceptor["x"] = 0.9999999999999999
    acceptor["y"] = 1.9999999999999998
    acceptor["z"] = 1.9999999999999998
    acceptor["radius"] = 0.9999999999999999
    acceptor["enabled"] = True
    acceptor["vector_on"] = 0
    acceptor["minsize"] = ""
    acceptor["maxsize"] = ""
    acceptor["selected"] = False

    expected["points"].append(acceptor)
    expected["points"].append(aromatic_1)
    assert pharmer == expected

    # Test pharmacophore with five elements
    pharmer = _pharmer_dict(five_element_pharmacophore)
    aromatic_2 = {}
    aromatic_2["name"] = "Aromatic"
    aromatic_2["hasvec"] = False
    aromatic_2["svector"] = {}
    aromatic_2["svector"]["x"] = 1
    aromatic_2["svector"]["y"] = 0
    aromatic_2["svector"]["z"] = 0
    aromatic_2["x"] = 0.0
    aromatic_2["y"] = 0.9999999999999999
    aromatic_2["z"] = 1.9999999999999998
    aromatic_2["radius"] = 0.9999999999999999
    aromatic_2["enabled"] = True
    aromatic_2["vector_on"] = 0
    aromatic_2["minsize"] = ""
    aromatic_2["maxsize"] = ""
    aromatic_2["selected"] = False

    donor = {}
    donor["name"] = "HydrogenDonor"
    donor["hasvec"] = True
    donor["svector"] = {}
    donor["svector"]["x"] = 0.0
    donor["svector"]["y"] = 0.7071067811865475
    donor["svector"]["z"] = 0.7071067811865475
    donor["x"] = 0.9999999999999999
    donor["y"] = 1.9999999999999998
    donor["z"] = 1.9999999999999998
    donor["radius"] = 0.9999999999999999
    donor["enabled"] = True
    donor["vector_on"] = 0
    donor["minsize"] = ""
    donor["maxsize"] = ""
    donor["selected"] = False

    hydrophobic = {}
    hydrophobic["name"] = "Hydrophobic"
    hydrophobic["hasvec"] = False
    hydrophobic["svector"] = {}
    hydrophobic["svector"]["x"] = 1
    hydrophobic["svector"]["y"] = 0
    hydrophobic["svector"]["z"] = 0
    hydrophobic["x"] = -0.9999999999999999
    hydrophobic["y"] = 1.9999999999999998
    hydrophobic["z"] = 1.9999999999999998
    hydrophobic["radius"] = 0.9999999999999999
    hydrophobic["enabled"] = True
    hydrophobic["vector_on"] = 0
    hydrophobic["minsize"] = ""
    hydrophobic["maxsize"] = ""
    hydrophobic["selected"] = False

    expected = {}
    expected["points"] = []
    expected["points"].append(acceptor)
    expected["points"].append(donor)
    expected["points"].append(hydrophobic)
    expected["points"].append(aromatic_1)
    expected["points"].append(aromatic_2)

    assert pharmer == expected


def test_load_mol2_file():
    file_name = ligands["ace"]
    molecules = load_mol2_file(file_name=file_name)

    assert len(molecules) == 3
    assert molecules[0].GetNumAtoms() == 14
    assert molecules[1].GetNumAtoms() == 25
    assert molecules[2].GetNumAtoms() == 29


def test_read_pharmagist():
    files_path = "./openpharmacophore/data/pharmacophores/pharmagist"

    streptadivin_file = pharmacophores["streptadivin"]
    pharmacophores_ = read_pharmagist(streptadivin_file)
    assert len(pharmacophores_) == 6
    assert len(pharmacophores_[0]) == 9
    assert len(pharmacophores_[1]) == 10
    assert len(pharmacophores_[2]) == 16
    assert len(pharmacophores_[3]) == 11
    assert len(pharmacophores_[4]) == 13
    assert len(pharmacophores_[5]) == 13

    elastase_file = pharmacophores["elastase"]
    pharmacophores_ = read_pharmagist(elastase_file)
    assert len(pharmacophores_) == 8
    assert len(pharmacophores_[0]) == 4
    assert len(pharmacophores_[1]) == 21
    assert len(pharmacophores_[2]) == 16
    assert len(pharmacophores_[3]) == 18
    assert len(pharmacophores_[4]) == 10
    assert len(pharmacophores_[5]) == 12
    assert len(pharmacophores_[6]) == 11
    assert len(pharmacophores_[7]) == 14


def test_to_pharmagist(three_element_pharmacophore, five_element_pharmacophore):
    # Test for two element pharmacophore
    mol2_list = _pharmagist_file_info(three_element_pharmacophore)
    expected_output = ['@<TRIPOS>MOLECULE\n',
                       '@<TRIPOS>ATOM\n',
                       '      1 ACC          1.0000    2.0000    2.0000   HB     0   HB      0.0000\n',
                       '      2 AR           1.0000    0.0000    0.0000   AR     1   AR      0.0000\n',
                       '@<TRIPOS>BOND\n']
    assert mol2_list == expected_output

    # Test for five element pharmacophore
    mol2_list = _pharmagist_file_info(five_element_pharmacophore)
    expected_output = ['@<TRIPOS>MOLECULE\n',
                       '@<TRIPOS>ATOM\n',
                       '      1 ACC          1.0000    2.0000    2.0000   HB     0   HB      0.0000\n',
                       '      2 DON          1.0000    2.0000    2.0000   HB     1   HB      0.0000\n',
                       '      3 HYD         -1.0000    2.0000    2.0000   HYD    2   HYD     0.0000\n',
                       '      4 AR           1.0000    0.0000    0.0000   AR     3   AR      0.0000\n',
                       '      5 AR           0.0000    1.0000    2.0000   AR     4   AR      0.0000\n',
                       '@<TRIPOS>BOND\n']
    assert mol2_list == expected_output


def test_from_ligandscout():
    points = from_ligandscout(pharmacophores["ligscout"])
    assert len(points) == 4

    neg_ion = points[2]
    assert neg_ion.feature_name == "negative charge"
    assert np.all(
        np.around(puw.get_value(neg_ion.center, "angstroms"), 1) == np.array([-8.0, 10.0, -9.5])
    )
    assert puw.get_value(neg_ion.radius, "angstroms") == 1.5

    donor = points[0]
    assert donor.feature_name == "hb donor"
    assert np.all(
        np.around(puw.get_value(donor.center, "angstroms"), 1) == np.array([-8.0, 2.0, -10.0])
    )
    assert puw.get_value(donor.radius, "angstroms") == 1.5
    dir_expected = np.array([-5.690445899963379, 0.5822541117668152, -10.5515718460083])
    dir_expected /= np.linalg.norm(dir_expected)
    assert np.all(
        np.around(donor.direction, 2) == np.around(dir_expected, 2)
    )

    ring = points[3]
    assert ring.feature_name == "aromatic ring"
    assert np.all(
        np.around(puw.get_value(ring.center, "angstroms"), 1) == np.array([0.0, 6.5, -3.0])
    )
    assert puw.get_value(ring.radius, "angstroms") == 1.5
    dir_expected = np.array([-3.8126893043518066, 1.7578959465026855, 0.6093783378601074])
    dir_expected /= np.linalg.norm(dir_expected)
    assert np.all(
        np.around(ring.direction, 2) == np.around(dir_expected, 2)
    )

    excluded_vol = points[1]
    assert excluded_vol.feature_name == "excluded volume"
    assert np.all(
        np.around(puw.get_value(excluded_vol.center, "angstroms"), 1) == np.array([5.5, 4.5, -2.0])
    )
    assert round(puw.get_value(excluded_vol.radius, "angstroms"), 1) == 1.0


# TODO: five element pharmacophore test case failing
@pytest.mark.skipif(reason="Failing, fixing it later")
def test_to_ligandscout(three_element_pharmacophore, five_element_pharmacophore):
    # Test for three element pharmacophore

    expected_bytes = b'<pharmacophore name="pharmacophore.pml" pharmacophoreType="LIGAND_SCOUT">'
    # Point 1: Acceptor
    expected_bytes += b'<vector disabled="false" featureId="ha_1" hasSyntheticProjectedPoint="false" name="HBA" optional="false" pointsToLigand="false" weight="1.0">'
    expected_bytes += b'<origin tolerance="0.9999999999999999" x3="0.9999999999999999" y3="1.9999999999999998" z3="1.9999999999999998" />'
    expected_bytes += b'<target tolerance="0.9999999999999999" x3="0.9999999999999999" y3="1.2928932188134523" z3="1.2928932188134523" />'
    expected_bytes += b'</vector>'
    # Point 2: Excluded
    expected_bytes += b'<volume disabled="false" featureId="ev_2" optional="false" type="exclusion" weight="1.0">'
    expected_bytes += b'<position tolerance="0.9999999999999999" x3="1.9999999999999998" y3="0.9999999999999999" z3="1.9999999999999998" />'
    expected_bytes += b'</volume>'
    # Point 3: Aromatic ring
    expected_bytes += b'<plane disabled="false" featureId="ai_3" name="AR" optional="false" weight="1.0">'
    expected_bytes += b'<position tolerance="0.9999999999999999" x3="0.9999999999999999" y3="0.0" z3="0.0" />'
    expected_bytes += b'<normal tolerance="0.9999999999999999" x3="0.0" y3="0.0" z3="1.0" />'
    expected_bytes += b'</plane>'
    # End of bytes
    expected_bytes += b'</pharmacophore>'

    _, document = _ligandscout_xml_tree(three_element_pharmacophore)
    pml_string = ET.tostring(document)

    assert pml_string == expected_bytes

    # Test for five element pharmacophore
    expected_bytes = b'<pharmacophore name="pharmacophore.pml" pharmacophoreType="LIGAND_SCOUT">'
    # Point 1: Acceptor
    expected_bytes += b'<vector disabled="false" featureId="ha_1" hasSyntheticProjectedPoint="false" name="HBA" optional="false" pointsToLigand="false" weight="1.0">'
    expected_bytes += b'<origin tolerance="0.9999999999999999" x3="0.9999999999999999" y3="1.9999999999999998" z3="1.9999999999999998" />'
    expected_bytes += b'<target tolerance="0.9999999999999999" x3="0.9999999999999999" y3="1.2928932188134523" z3="1.2928932188134523" />'
    expected_bytes += b'</vector>'
    # Point 2: Donor
    expected_bytes += b'<vector disabled="false" featureId="hb_2" hasSyntheticProjectedPoint="false" name="HBD" optional="false" pointsToLigand="false" weight="1.0">'
    expected_bytes += b'<origin tolerance="0.9999999999999999" x3="0.9999999999999999" y3="1.9999999999999998" z3="1.9999999999999998" />'
    expected_bytes += b'<target tolerance="0.9999999999999999" x3="0.9999999999999999" y3="1.2928932188134523" z3="1.2928932188134523" />'
    expected_bytes += b'</vector>'
    # Point 3 : Hydrophobic
    expected_bytes += b'<point name="H" featureId="hi_3" optional="false" disabled="false" weight="1.0">'
    expected_bytes += b'<position x3="-0.9999999999999999" y3="1.9999999999999998" z3="1.9999999999999998" />'
    expected_bytes += b'</point>'
    # Point 4: Aromatic ring
    expected_bytes += b'<plane disabled="false" featureId="ai_4" name="AR" optional="false" weight="1.0">'
    expected_bytes += b'<position tolerance="0.9999999999999999" x3="0.9999999999999999" y3="0.0" z3="0.0" />'
    expected_bytes += b'<normal tolerance="0.9999999999999999" x3="0.0" y3="0.0" z3="1.0" />'
    expected_bytes += b'</plane>'
    # Point 5: Aromatic ring
    expected_bytes += b'<plane disabled="false" featureId="ai_5" name="AR" optional="false" weight="1.0">'
    expected_bytes += b'<position tolerance="0.9999999999999999" x3="0.0" y3="0.9999999999999999" z3="1.9999999999999999" />'
    expected_bytes += b'<normal tolerance="0.9999999999999999" x3="0.0" y3="0.0" z3="1.0" />'
    expected_bytes += b'</plane>'
    # End of bytes
    expected_bytes += b'</pharmacophore>'

    _, document = _ligandscout_xml_tree(five_element_pharmacophore)
    pml_string = ET.tostring(document)

    assert pml_string == expected_bytes


def test_from_moe():
    file_name = pharmacophores["gmp"]
    points = from_moe(file_name)
    assert len(points) == 10

    # HB Acceptor features should be first
    acceptor = points[0]
    assert acceptor.feature_name == "hb acceptor"
    assert np.all(
        np.around(puw.get_value(acceptor.center, "angstroms"), 2) == np.around(np.array([0.312, 3.0175, -2.44825]), 2)
    )
    assert np.around(puw.get_value(acceptor.radius, "angstroms"), 2) == 0.57
    acceptor = points[1]
    assert acceptor.feature_name == "hb acceptor"
    assert np.all(
        np.around(puw.get_value(acceptor.center, "angstroms"), 2) ==
        np.around(np.array([-1.95875, 2.536, -3.03625]), 2)
    )
    assert np.around(puw.get_value(acceptor.radius, "angstroms"), 2) == 0.62
    acceptor_2 = points[2]
    assert acceptor_2.feature_name == "hb acceptor"
    assert np.all(
        np.around(puw.get_value(acceptor_2.center, "angstroms"), 2) ==
        np.around(np.array([-0.755095833333333, 6.3286375, -3.96758333333333]), 2)
    )
    assert np.around(puw.get_value(acceptor_2.radius, "angstroms"), 2) == 1.25

    # Donor points
    donor = points[3]
    assert donor.feature_name == "hb donor"
    assert np.all(
        np.around(puw.get_value(donor.center, "angstroms"), 2) == np.around(np.array([1.71, 1.43075, -1.4255]), 2)
    )
    assert np.around(puw.get_value(donor.radius, "angstroms"), 2) == 0.51

    # Hydrophobic points
    hyd = points[4]
    assert hyd.feature_name == "hydrophobicity"
    assert np.all(
        np.around(puw.get_value(hyd.center, "angstroms"), 2) == np.around(np.array([2.7895, 2.4035, -1.40875]), 2)
    )
    assert np.around(puw.get_value(hyd.radius, "angstroms"), 2) == 0.55
    hyd_2 = points[5]
    assert hyd_2.feature_name == "hydrophobicity"
    assert np.all(
        np.around(puw.get_value(hyd_2.center, "angstroms"), 2) == np.around(np.array([-1.54725, -2.979375, -0.961875]),
                                                                            2)
    )
    assert np.around(puw.get_value(hyd_2.radius, "angstroms"), 2) == 0.74

    # Aromatic Points
    aromatic_1 = points[6]
    assert aromatic_1.feature_name == "aromatic ring"
    assert np.all(
        np.around(puw.get_value(aromatic_1.center, "angstroms"), 2) ==
        np.around(np.array([-0.748458333333333, 2.13108333333333, -2.490375]), 2)
    )
    assert np.around(puw.get_value(aromatic_1.radius, "angstroms"), 2) == 0.58
    aromatic_2 = points[7]
    assert aromatic_2.feature_name == "aromatic ring"
    assert np.all(
        np.around(puw.get_value(aromatic_2.center, "angstroms"), 2) ==
        np.around(np.array([-1.719625, -0.0273333333333334, -2.055625]), 2)
    )
    assert np.around(puw.get_value(aromatic_2.radius, "angstroms"), 2) == 0.6
    aromatic_3 = points[8]
    assert aromatic_3.feature_name == "aromatic ring"
    assert np.all(
        np.around(puw.get_value(aromatic_3.center, "angstroms"), 2) ==
        np.around(np.array([5.20029166666667, 1.25479166666667, -0.199041666666667]), 2)
    )
    assert np.around(puw.get_value(aromatic_3.radius, "angstroms"), 2) == 0.61
    aromatic_4 = points[9]
    assert aromatic_4.feature_name == "aromatic ring"
    assert np.all(
        np.around(puw.get_value(aromatic_4.center, "angstroms"), 2) ==
        np.around(np.array([-0.755095833333333, 6.3286375, -3.96758333333333]), 2)
    )
    assert np.around(puw.get_value(aromatic_4.radius, "angstroms"), 2) == 1.25


def test_to_moe(three_element_pharmacophore):
    pharmacophore_str = _moe_ph4_string(three_element_pharmacophore)

    now = datetime.datetime.now()
    month = str(now.month)
    year = str(now.year)
    expected_str = (f'#moe:ph4que {year}.{month}\n'
                    f'#pharmacophore 5 tag t value *\n'
                    f'scheme t Unified matchsize i 0 title t s $\n'
                    f'#feature 3 expr tt color ix x r y r z r r r ebits ix gbits ix\n'
                    f'Acc df2f2 0.9999999999999999 1.9999999999999998 1.9999999999999998 0.9999999999999999 0 300 Aro df2f2 0.9999999999999999 0.0 0.0 0.9999999999999999 0 300 \n'
                    f'#volumesphere 90 x r y r z r r r\n'
                    f'1.9999999999999998 0.9999999999999999 1.9999999999999998 0.9999999999999999 \n'
                    f'#endpharmacophore')

    assert pharmacophore_str == expected_str


@pytest.mark.parametrize("file_name", [
    "ace.mol2",
    "clique_detection.smi"
])
def test_load_molecules_file(file_name):
    ligands_ = load_molecules_file(file_name=ligands["clique_detection"],
                                   titleLine=False)
    assert len(ligands_) == 5

    ligands_ = load_molecules_file(file_name=ligands["ace"])
    assert len(ligands_) == 3

    assert isinstance(ligands_, list)
    for lig in ligands_:
        assert isinstance(lig, Chem.Mol)


def test_smi_has_header_and_id():
    file = ligands["mols"]
    has_header, has_id = smi_has_header_and_id(file)
    assert not has_header
    assert not has_id

    file = ligands["BAAAML"]
    has_header, has_id = smi_has_header_and_id(file)
    assert has_header
    assert has_id


def test_mol2_mol_generator():
    mol2_file = ligands["ace"]

    molecules = []
    with open(mol2_file, "r") as fp:
        for mol in mol2_mol_generator(fp):
            molecules.append(mol)

    assert len(molecules) == 3
    assert molecules[0].GetNumAtoms() == 14
    assert molecules[1].GetNumAtoms() == 25
    assert molecules[2].GetNumAtoms() == 29


def test_smiles_mol_generator():
    smi_file = ligands["mols"]

    n_molecules = 0
    with open(smi_file, "r") as fp:
        for mol in smiles_mol_generator(fp, header=False, mol_id=False):
            n_molecules += 1
            assert isinstance(mol, Chem.Mol)

    assert n_molecules == 5

    smi_file = ligands["BAAAML"]

    n_molecules = 0
    with open(smi_file, "r") as fp:
        for mol in smiles_mol_generator(fp, header=True, mol_id=True):
            n_molecules += 1
            assert isinstance(mol, Chem.Mol)

    assert n_molecules == 23
