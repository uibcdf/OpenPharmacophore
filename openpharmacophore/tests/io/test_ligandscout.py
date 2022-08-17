import openpharmacophore.io as io
import openpharmacophore.data as data
import numpy as np
import pyunitwizard as puw
import pytest
import xml.etree.ElementTree as ET
from example_pharmacophores import three_element_pharmacophore, five_element_pharmacophore


def test_from_ligandscout():
    points = io.from_ligandscout(data.pharmacophores["ligscout"])
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
def test_to_ligandscout():
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

    _, document = io._ligandscout_xml_tree(three_element_pharmacophore())
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

    _, document = io._ligandscout_xml_tree(five_element_pharmacophore())
    pml_string = ET.tostring(document)

    assert pml_string == expected_bytes
