from openpharmacophore import PharmacophoricPoint
from collections import namedtuple
import xml.etree.ElementTree as ET
import pyunitwizard as puw
import numpy as np

def from_ligandscout(file_name): 
    """ Loads a pharmacophore from a ligandscout pml file.

        Parameters
        ----------
        pharmacophore_file : str
            Name of the file containing the pharmacophore.

        Returns
        -------
        points : list of openpharmacophore.PharmacophoricPoints
            A list of pharmacophoric elements.
        
    """
    tree = ET.parse(file_name)
    pharmacophore = tree.getroot()

    ligandscout_to_oph = { # dictionary to map ligandscout features to openpharmacophore elements
        "PI": "positive charge",
        "NI": "negative charge",
        "HBD": "hb donor", 
        "HBA": "hb acceptor",
        "H": "hydrophobicity",
        "AR": "aromatic ring",
        "exclusion": "excluded volume"
    }

    points = []
    for element in pharmacophore:

        if element.tag == "point":
            feat_name = element.attrib["name"]
            for position in element:
                coords = position.attrib
                x = float(coords["x3"])
                y = float(coords["y3"])
                z = float(coords["z3"])
                radius = float(coords["tolerance"])
            point = PharmacophoricPoint(feat_type=ligandscout_to_oph[feat_name],
                        center=puw.quantity([x, y, z], "angstroms"),
                        radius=puw.quantity(radius, "angstroms"))
            points.append(point)

        elif element.tag == "vector":
            feat_name = element.attrib["name"]
            for position in element:
                if position.tag == "origin":
                    coords_1 = position.attrib
                    x_1 = float(coords_1["x3"])
                    y_1 = float(coords_1["y3"])
                    z_1 = float(coords_1["z3"])
                    radius_1 = float(coords_1["tolerance"])
                elif position.tag == "target":
                    coords_2 = position.attrib
                    x_2 = float(coords_2["x3"])
                    y_2 = float(coords_2["y3"])
                    z_2 = float(coords_2["z3"])
                    radius_2 = float(coords_2["tolerance"])  
            point = PharmacophoricPoint(feat_type=ligandscout_to_oph[feat_name],
                        center=puw.quantity([x_1, y_1, z_1], "angstroms"),
                        radius=puw.quantity(radius_1, "angstroms"),
                        direction=[x_2, y_2, z_2])
            points.append(point)

        elif element.tag == "plane":
            feat_name = element.attrib["name"]
            for position in element:
                if position.tag == "position":
                    coords_1 = position.attrib
                    x_1 = float(coords_1["x3"])
                    y_1 = float(coords_1["y3"])
                    z_1 = float(coords_1["z3"])
                    radius_1 = float(coords_1["tolerance"])
                    center = np.array([x_1, y_1, z_1])
                elif position.tag == "normal":
                    coords_2 = position.attrib
                    x_2 = float(coords_2["x3"])
                    y_2 = float(coords_2["y3"])
                    z_2 = float(coords_2["z3"])
                    radius_2 = float(coords_2["tolerance"])
            point = PharmacophoricPoint(
                        feat_type=ligandscout_to_oph[feat_name],
                        center=puw.quantity(center, "angstroms"),
                        radius=puw.quantity(radius_1, "angstroms"),
                        direction=[x_2, y_2, z_2])
            points.append(point)

        elif element.tag == "volume":
            feat_name = element.attrib["type"]
            for position in element:
                coords = position.attrib
                x = float(coords["x3"])
                y = float(coords["y3"])
                z = float(coords["z3"])
                radius = float(coords["tolerance"])
            point = PharmacophoricPoint(
                    feat_type=ligandscout_to_oph[feat_name],
                    center=puw.quantity([x, y, z], "angstroms"),
                    radius=puw.quantity(radius, "angstroms"))
            points.append(point)

    return points

def to_ligandscout(pharmacophore, file_name):
    """ Save a pharmacophore as a ligandscout file (pml file)

        Parameters
        ----------
        pharmacophore : openpharmacophore.Pharmacophore
            Pharmacophore object that will be saved to a file. Can be a Pharmacophore, 
            LigandBasedPharmacophore or StructureBasedPharmaophore.

        file_name : str
            Name of the file that will contain the pharmacophore.

        Note
        ----
        Nothing is returned. A new file is written.
    """
    tree, _ = _ligandscout_xml_tree(pharmacophore)
    tree.write(file_name, encoding="UTF-8", xml_declaration=True)

def _ligandscout_xml_tree(pharmacophore):
    """ Get an xml element tree necesary to create a ligandscout pharmacophore.

        Parameters
        ----------
        pharmacophore : openpharmacophore.Pharmacophore
            Pharmacophore object that will be saved to a file.

        Returns
        -------
        tree : xml.etree.ElementTree
            The element tree.

    """
    Feature = namedtuple("Feature", ["name", "id"])
    feature_mapper = { # dictionary to map openpharmacophore features to ligandscout
        "aromatic ring": Feature("AR", "ai_"),
        "hydrophobicity": Feature("H", "hi_"),
        "hb acceptor": Feature("HBA", "ha_"),
        "hb donor": Feature("HBD", "hd_"),
        "excluded volume": Feature("exclusion", "ev_"),
        "positive charge": Feature("PI", "pi_"),
        "negative charge": Feature("NI", "ni_"),
    }

    tree = ET.ElementTree("tree")
    document = ET.Element("pharmacophore")
    document.set("name", "pharmacophore.pml")
    document.set("pharmacophoreType", "LIGAND_SCOUT")

    for i, element in enumerate(pharmacophore.elements):
        try:
            feat_name = feature_mapper[element.feature_name].name
        except: # skip features not supported by ligandscout
            continue
        coords = puw.get_value(element.center, to_unit="angstroms")
        x = str(coords[0])
        y = str(coords[1])
        z = str(coords[2])
        radius = str(puw.get_value(element.radius, to_unit="angstroms"))
        feat_id =  feature_mapper[element.feature_name].id + str(i + 1)

        is_point = (feat_name == "PI" or feat_name == "NI" or feat_name == "H" 
                    or feat_name == "HBD" or feat_name == "HBA")
        is_vector = element.has_direction and (feat_name == "HBD" or feat_name == "HBA")

        if is_vector:
            direction = coords - element.direction
            dir_x = str(direction[0])
            dir_y = str(direction[1])
            dir_z = str(direction[2])
            vector = ET.SubElement(document, "vector")
            # Set vector attributes
            vector.set("name", feat_name)
            vector.set("featureId", feat_id)
            vector.set("pointsToLigand", "false")
            vector.set("hasSyntheticProjectedPoint", "false")
            vector.set("optional", "false")
            vector.set("disabled", "false")
            vector.set("weight", "1.0")
            # Add origin tag
            origin = ET.SubElement(vector, "origin")
            origin.set("x3", x)
            origin.set("y3", y)
            origin.set("z3", z)
            origin.set("tolerance", radius)
            # Add target tag
            target = ET.SubElement(vector, "target")
            target.set("x3", dir_x)
            target.set("y3", dir_y)
            target.set("z3", dir_z)
            target.set("tolerance", radius)
        elif is_point:
            point = ET.SubElement(document, "point")
            # Set point attributes
            point.set("name", feat_name)
            point.set("featureId", feat_id)
            point.set("optional", "false")
            point.set("disabled", "false")
            point.set("weight", "1.0")
            # Add position tag
            position = ET.SubElement(point, "position")
            position.set("x3", x)
            position.set("y3", y)
            position.set("z3", z)
            position.set("tolerance", radius)
        elif feat_name == "AR":
            direction = element.direction
            dir_x = str(direction[0])
            dir_y = str(direction[1])
            dir_z = str(direction[2])
            plane = ET.SubElement(document, "plane")
            # Set plane attributes
            plane.set("name", feat_name)
            plane.set("featureId", feat_id)
            plane.set("optional", "false")
            plane.set("disabled", "false")
            plane.set("weight", "1.0")
            # Add position tag
            position = ET.SubElement(plane, "position")
            position.set("x3", x)
            position.set("y3", y)
            position.set("z3", z)
            position.set("tolerance", radius)
            # Add normal tag
            normal = ET.SubElement(plane, "normal")
            normal.set("x3", dir_x)
            normal.set("y3", dir_y)
            normal.set("z3", dir_z)
            normal.set("tolerance", radius)
        elif feat_name == "exclusion":
            volume = ET.SubElement(document, "volume")
            # Set volume attributes
            volume.set("type", "exclusion")
            volume.set("featureId", feat_id)
            volume.set("optional", "false")
            volume.set("disabled", "false")
            volume.set("weight", "1.0")
            # Add position tag
            position = ET.SubElement(volume, "position")
            position.set("x3", x)
            position.set("y3", y)
            position.set("z3", z)
            position.set("tolerance", radius)

    tree._setroot(document)

    return tree, document
