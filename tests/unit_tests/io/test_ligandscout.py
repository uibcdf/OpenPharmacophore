from openpharmacophore.io.ligandscout import read_ligandscout, ligandscout_xml_tree
import openpharmacophore.data as data
from example_pharmacophores import three_element_pharmacophore, five_element_pharmacophore
import numpy as np
import pyunitwizard as puw
import math
import xml.etree.ElementTree as ET


def test_from_ligandscout():
    points = read_ligandscout(data.pharmacophores["ligscout"])
    assert len(points) == 4

    neg_ion = points[0]
    assert neg_ion.feature_name == "negative charge"
    assert np.all(
        np.around(puw.get_value(neg_ion.center, "angstroms"), 1) == np.array([-8.0, 10.0, -9.5])
    )
    assert puw.get_value(neg_ion.radius, "angstroms") == 1.5

    donor = points[1]
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

    ring = points[2]
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

    excluded_vol = points[3]
    assert excluded_vol.feature_name == "excluded volume"
    assert np.all(
        np.around(puw.get_value(excluded_vol.center, "angstroms"), 1) == np.array([5.5, 4.5, -2.0])
    )
    assert round(puw.get_value(excluded_vol.radius, "angstroms"), 1) == 1.0


def test_to_ligandscout_with_tree_element_pharmacophore():
    # Test for three element pharmacophore

    expected_bytes = b'<pharmacophore name="pharmacophore.pml" pharmacophoreType="LIGAND_SCOUT">'
    # Point 1: Acceptor
    expected_bytes += b'<vector disabled="false" featureId="ha_1" hasSyntheticProjectedPoint="false"'
    expected_bytes += b' name="HBA" optional="false" pointsToLigand="false" weight="1.0">'
    expected_bytes += b'<origin tolerance="1.000" x3="1.000" y3="2.000" z3="2.000" />'
    expected_bytes += b'<target tolerance="1.000" x3="1.000" y3="1.293" z3="1.293" />'
    expected_bytes += b'</vector>'
    # Point 2: Excluded
    expected_bytes += b'<volume disabled="false" featureId="ev_2" optional="false" type="exclusion" weight="1.0">'
    expected_bytes += b'<position tolerance="1.000" x3="2.000" y3="1.000" z3="2.000" />'
    expected_bytes += b'</volume>'
    # Point 3: Aromatic ring
    expected_bytes += b'<plane disabled="false" featureId="ai_3" name="AR" optional="false" weight="1.0">'
    expected_bytes += b'<position tolerance="1.000" x3="1.000" y3="0.000" z3="0.000" />'
    expected_bytes += b'<normal tolerance="1.000" x3="0.000" y3="0.000" z3="1.000" />'
    expected_bytes += b'</plane>'
    # End of bytes
    expected_bytes += b'</pharmacophore>'

    tree = ligandscout_xml_tree(three_element_pharmacophore())
    pml_string = ET.tostring(tree.getroot())

    assert pml_string == expected_bytes


def test_ligand_scout_xml_tree():

    tree = ligandscout_xml_tree(five_element_pharmacophore())
    document = tree.getroot()

    assert len(document) == 5
    assert document.tag == "pharmacophore"

    vectors = document.findall("vector")
    assert len(vectors) == 2
    assert vectors[0].get("name") == "HBA"
    assert vectors[0].get("featureId") == "ha_1"
    origin = vectors[0].findall("origin")
    assert len(origin) == 1
    assert math.ceil(float(origin[0].get("x3"))) == 1.0
    assert math.ceil(float(origin[0].get("y3"))) == 2.0
    assert math.ceil(float(origin[0].get("z3"))) == 2.0
    target = vectors[0].findall("target")
    assert len(target) == 1
    assert math.ceil(float(target[0].get("x3"))) == 1.0
    assert round(float(target[0].get("y3")), 1) == 1.3
    assert round(float(target[0].get("z3")), 1) == 1.3

    assert vectors[1].get("name") == "HBD"
    assert vectors[1].get("featureId") == "hd_2"

    hydrophobic = document.findall("point")
    assert len(hydrophobic) == 1
    assert hydrophobic[0].get("name") == "H"
    assert hydrophobic[0].get("featureId") == "hi_3"
    position = hydrophobic[0].findall("position")
    assert len(position) == 1
    assert math.floor(float(position[0].get("x3"))) == -1.0
    assert math.ceil(float(position[0].get("y3"))) == 2.0
    assert math.ceil(float(position[0].get("z3"))) == 2.0

    rings = document.findall("plane")
    assert len(rings) == 2
    ring_1 = rings[0]
    ring_2 = rings[1]
    assert ring_1.get("name") == "AR"
    assert ring_1.get("featureId") == "ai_4"
    assert ring_2.get("name") == "AR"
    assert ring_2.get("featureId") == "ai_5"

    position = ring_1.findall("position")
    assert len(position) == 1
    assert math.ceil(float(position[0].get("x3"))) == 1.0
    assert math.ceil(float(position[0].get("y3"))) == 0.0
    assert math.ceil(float(position[0].get("z3"))) == 0.0
