import datetime
import math
import xml.etree.ElementTree as ET
from openpharmacophore import pharmacophore_writer


class TestWriteMol2:

    def test_pad_coordinate_with_zeros_positive_number(self):
        assert pharmacophore_writer._pad_coordinate_with_zeros(25) == "25.000"
        assert pharmacophore_writer._pad_coordinate_with_zeros(2.5) == "2.5000"
        assert pharmacophore_writer._pad_coordinate_with_zeros(0.25) == "0.2500"

    def test_pad_coordinate_with_zeros_negative_number(self):
        assert pharmacophore_writer._pad_coordinate_with_zeros(-25) == "-25.000"
        assert pharmacophore_writer._pad_coordinate_with_zeros(-2.5) == "-2.5000"
        assert pharmacophore_writer._pad_coordinate_with_zeros(-0.25) == "-0.2500"

    def test_pad_number_with_zeros_integer(self):
        assert pharmacophore_writer._pad_coordinate_with_zeros(1) == "1.0000"
        assert pharmacophore_writer._pad_coordinate_with_zeros(1.0) == "1.0000"

    def test_mol2_pharmacophores(
            self,
            two_element_pharmacophore,
            five_element_pharmacophore
    ):
        # Test for two element pharmacophore
        ph_1 = two_element_pharmacophore
        ph_1.score = 0.8452
        ph_1.ref_mol = 2
        ph_1.ref_struct = 4
        ph_1.props = {"min_actives": 5}

        mol2_list = pharmacophore_writer._mol2_pharmacophores([ph_1], save_props=True)
        expected_output_1 = [
            '@<TRIPOS>PHARMACOPHORE\n',
            '@<TRIPOS>POINTS\n',
            '      1 AR           1.0000    0.0000    0.0000   AR     0   AR      1.0000\n',
            '      2 ACC          1.0000    2.0000    2.0000   HB     1   HB      1.0000\n',
            '@<TRIPOS>PROPERTIES\n',
            '       min_actives        5\n',
            '       score         0.8452\n',
            '       ref_mol            2\n',
            '       ref_struct         4\n'
        ]
        assert mol2_list == expected_output_1

        # Test for five element pharmacophore
        ph_2 = five_element_pharmacophore
        mol2_list = pharmacophore_writer._mol2_pharmacophores([ph_2], save_props=True)
        expected_output_2 = [
            '@<TRIPOS>PHARMACOPHORE\n',
            '@<TRIPOS>POINTS\n',
            '      1 ACC          1.0000    2.0000    2.0000   HB     0   HB      1.0000\n',
            '      2 DON          1.0000    2.0000    2.0000   HB     1   HB      1.0000\n',
            '      3 HYD         -1.0000    2.0000    2.0000   HYD    2   HYD     1.0000\n',
            '      4 AR           1.0000    0.0000    0.0000   AR     3   AR      1.0000\n',
            '      5 AR           0.0000    1.0000    2.0000   AR     4   AR      1.0000\n',
        ]
        assert mol2_list == expected_output_2

        # Test for multiple pharmacophores
        expected_output_3 = expected_output_1 + expected_output_2
        mol2_list = pharmacophore_writer._mol2_pharmacophores([ph_1, ph_2], save_props=True)
        assert mol2_list == expected_output_3


class TestWritePH4:

    def test_ph4_pharmacophore(self, three_element_pharmacophore):
        pharmacophore_str = pharmacophore_writer._ph4_pharmacophore(three_element_pharmacophore)

        now = datetime.datetime.now()
        month = str(now.month)
        year = str(now.year)
        expected_str = (f'#moe:ph4que {year}.{month}\n'
                        f'#pharmacophore 5 tag t value *\n'
                        f'scheme t Unified matchsize i 0 title t s $\n'
                        f'#feature 3 expr tt color ix x r y r z r r r ebits ix gbits ix\n'
                        f'Acc df2f2 1.000 2.000 2.000 1.000 0 300 Aro df2f2 1.000 0.000 0.000 1.000 0 300 \n'
                        f'#volumesphere 90 x r y r z r r r\n'
                        f'2.000 1.000 2.000 1.000 \n'
                        f'#endpharmacophore')

        assert pharmacophore_str == expected_str


class TestWriteJSON:

    @staticmethod
    def create_json_pharmacophoric_element(name, coords, radius, direction=None):

        has_direction = direction is not None
        if has_direction:
            direction = {
                "x": direction[0],
                "y": direction[1],
                "z": direction[2],
            }
        else:
            direction = {
                "x": 1.,
                "y": 0.,
                "z": 0.,
            }

        return {
            "name": name,
            "hasvec": has_direction,
            "svector": direction,
            "x": coords[0],
            "y": coords[1],
            "z": coords[2],
            "radius": radius,
            "enabled": True,
            "vector_on": 0,
            "minsize": "",
            "maxsize": "",
            "selected": False
        }

    def test_json_pharmacophoric_elements(self, three_element_pharmacophore):

        # Expected output from to_json
        expected = {"points": []}

        aromatic = self.create_json_pharmacophoric_element(
            name="Aromatic",
            coords=(1.0, 0.0, 0.0),
            radius=1.0,
            direction=(0.0, 0.0, 1.0)
        )
        acceptor = self.create_json_pharmacophoric_element(
            name="HydrogenAcceptor",
            coords=(1.0, 2.0, 2.0),
            radius=1.0,
            direction=(0.0, 1 / math.sqrt(2.0), 1 / math.sqrt(2.0))
        )
        excluded = self.create_json_pharmacophoric_element(
            name="ExclusionSphere",
            coords=(2.0, 1.0, 2.0),
            radius=1.0
        )

        expected["points"].append(acceptor)
        expected["points"].append(excluded)
        expected["points"].append(aromatic)

        assert pharmacophore_writer._json_pharmacophore(three_element_pharmacophore) == expected


class TestWriteLigandScout:

    def test_to_ligandscout_with_tree_element_pharmacophore(self, three_element_pharmacophore):
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

        tree = pharmacophore_writer._ligandscout_xml_tree(three_element_pharmacophore)
        pml_string = ET.tostring(tree.getroot())

        assert pml_string == expected_bytes

    def test_ligand_scout_xml_tree(self, five_element_pharmacophore):
        tree = pharmacophore_writer._ligandscout_xml_tree(five_element_pharmacophore)
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
