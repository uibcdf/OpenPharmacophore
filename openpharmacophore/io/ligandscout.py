from collections import namedtuple
import pyunitwizard as puw
import xml.etree.ElementTree as ET
from ..pharmacophore import PharmacophoricPoint

ligandscout_to_oph = {
    "PI": "positive charge",
    "NI": "negative charge",
    "HBD": "hb donor",
    "HBA": "hb acceptor",
    "H": "hydrophobicity",
    "AR": "aromatic ring",
    "exclusion": "excluded volume"
}


def pharmacophoric_point_from_point_or_volume(element):
    """ Creates a pharmacophoric point from a ligandscout point or volume

        Returns
        -------
        PharmacophoricPoint
            The point
    """
    if element.tag == "point":
        feat_name = element.attrib["name"]
    else:
        feat_name = element.attrib["type"]
    coords = element[0].attrib
    x = float(coords["x3"])
    y = float(coords["y3"])
    z = float(coords["z3"])
    radius = float(coords["tolerance"])

    return PharmacophoricPoint(ligandscout_to_oph[feat_name],
                               puw.quantity([x, y, z], "angstroms"),
                               puw.quantity(radius, "angstroms"))


def pharmacophoric_point_from_vector_or_plane(element):
    """ Creates a pharmacophoric point from a ligandscout plane or vector

        Returns
        -------
        PharmacophoricPoint
            The point
        """
    feat_name = element.attrib["name"]
    if element.tag == "vector":
        tags = ["origin", "target"]
    else:
        tags = ["position", "normal"]

    for position in element:
        if position.tag == tags[0]:
            coords_1 = position.attrib
            x_1 = float(coords_1["x3"])
            y_1 = float(coords_1["y3"])
            z_1 = float(coords_1["z3"])
            radius_1 = float(coords_1["tolerance"])
        elif position.tag == tags[1]:
            coords_2 = position.attrib
            x_2 = float(coords_2["x3"])
            y_2 = float(coords_2["y3"])
            z_2 = float(coords_2["z3"])
    return PharmacophoricPoint(ligandscout_to_oph[feat_name],
                               puw.quantity([x_1, y_1, z_1], "angstroms"),
                               puw.quantity(radius_1, "angstroms"),
                               direction=[x_2, y_2, z_2])


def read_ligandscout(file_name):
    """ Reads a ligandscout pharmacophore (pml) file and returns a list of pharmacophoric points.

        Parameters
        ----------
        file_name : str
            Name of the file containing the pharmacophore.

        Returns
        -------
        points : list[PharmacophoricPoint]
            A list of pharmacophoric pharmacophoric_points.
        
    """
    tree = ET.parse(file_name)
    pharmacophore = tree.getroot()

    points = []
    for element in pharmacophore:

        if element.tag == "point" or element.tag == "volume":
            point = pharmacophoric_point_from_point_or_volume(element)

        elif element.tag == "vector" or element.tag == "plane":
            point = pharmacophoric_point_from_vector_or_plane(element)

        else:
            continue

        points.append(point)

    return points


def set_coords_and_radius(tree_element, coords, radius,
                          feat_id):
    """ Set the coords and radius of a tree element corresponding
        to a volume, plane or point tag.
    """
    x, y, z = coords

    tree_element.set("featureId", feat_id)
    tree_element.set("optional", "false")
    tree_element.set("disabled", "false")
    tree_element.set("weight", "1.0")
    # Add position tag
    position = ET.SubElement(tree_element, "position")
    position.set("x3", x)
    position.set("y3", y)
    position.set("z3", z)
    position.set("tolerance", radius)


def ligandscout_xml_tree(pharmacophoric_points):
    """ Returns a xml element tree necessary to create a ligandscout pharmacophore file.

        Parameters
        ----------
        pharmacophoric_points : openpharmacophore.Pharmacophore
            The pharmacophore that will be saved to a file.

        Returns
        -------
        tree : xml.etree.ElementTree
            The element tree.

    """
    Feature = namedtuple("Feature", ["name", "id"])
    feature_mapper = {
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

    for ii, element in enumerate(pharmacophoric_points):
        try:
            feat_name = feature_mapper[element.feature_name].name
        except KeyError:  # skip features not supported by ligandscout
            continue
        coords = puw.get_value(element.center, to_unit="angstroms")
        x = f"{coords[0]:.3f}"
        y = f"{coords[1]:.3f}"
        z = f"{coords[2]:.3f}"
        radius_val = puw.get_value(element.radius, to_unit="angstroms")
        radius = f"{radius_val:.3f}"
        feat_id = feature_mapper[element.feature_name].id + str(ii + 1)

        is_point = (feat_name == "PI" or feat_name == "NI" or feat_name == "H"
                    or feat_name == "HBD" or feat_name == "HBA")
        is_vector = element.has_direction and (feat_name == "HBD" or feat_name == "HBA")

        if is_vector:
            direction = coords - element.direction
            dir_x = f"{direction[0]:.3f}"
            dir_y = f"{direction[1]:.3f}"
            dir_z = f"{direction[2]:.3f}"
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
            point.set("name", feat_name)
            set_coords_and_radius(point, (x, y, z),
                                  radius, feat_id)
        elif feat_name == "AR":
            # TODO: ligandscout expects aromatic rings to have direction.
            #  what should we do if an aromatic point doesn't have direction
            if not element.has_direction:
                continue
            direction = element.direction
            dir_x = f"{direction[0]:.3f}"
            dir_y = f"{direction[1]:.3f}"
            dir_z = f"{direction[2]:.3f}"
            plane = ET.SubElement(document, "plane")
            plane.set("name", feat_name)
            set_coords_and_radius(plane, (x, y, z),
                                  radius, feat_id)
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
            set_coords_and_radius(volume, (x, y, z),
                                  radius, feat_id)

    tree._setroot(document)

    return tree
