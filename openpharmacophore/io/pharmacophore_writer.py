from collections import namedtuple
from copy import deepcopy
import datetime
import json
import xml.etree.ElementTree as ET
import pyunitwizard as puw
import numpy as np


# MOL2 Files

def _pad_coordinate_with_zeros(coord):
    """ Pad a number with zeros to the right.

        Parameters
        ----------
        coord: float
            A number

        Returns
        -------
        str
            A string with the coordinate padded with zeros to the right
    """
    # The final string will contain 6 characters in total including the point or 7
    # if the number is negative
    total_characters = 6
    coord_float = float(coord)

    if coord_float.is_integer():
        coord_str = str(coord) + "."
    else:
        coord_str = str(coord)

    if coord_float < 0:
        # If the coordinate is negative the sign will add a character
        total_characters += 1
    return coord_str.ljust(total_characters, "0")


def _get_props_dict(pharmacophore):
    """
        Parameters
        ----------
        pharmacophore: Pharmacophore

        Returns
        -------
        dict[str, str]
    """
    props = deepcopy(pharmacophore.props)
    if pharmacophore.score is not None:
        props["score"] = str(pharmacophore.score)
    if pharmacophore.ref_mol is not None:
        props["ref_mol"] = str(pharmacophore.ref_mol)
    if pharmacophore.ref_struct is not None:
        props["ref_struct"] = str(pharmacophore.ref_struct)
    return props


def _mol2_pharmacophores(pharmacophores, save_props):
    """ Get necessary info to create a mol2 file to store pharmacophores.

        Parameters
        ----------
        pharmacophores : list[Pharmacophore]
            Nested list of pharmacophoric points were each entry represents a pharmacophore

        save_props : bool
            Whether to save the properties stored in the pharmacophore.

        Returns
        -------
        lines : list[str]
            List with the contents of a mol2 file.
    """
    pharmagist_element_name = {
        "aromatic ring": "AR ",
        "hydrophobicity": "HYD",
        "hb acceptor": "ACC",
        "hb donor": "DON",
        "positive charge": "CAT",
        "negative charge": "ANI",
    }

    pharmagist_element_specs = {
        "aromatic ring": "AR ",
        "hydrophobicity": "HYD",
        "hb acceptor": "HB ",
        "hb donor": "HB ",
        "positive charge": "CHG",
        "negative charge": "CHG",
    }

    lines = []
    for pharmacophore in pharmacophores:
        lines.extend(["@<TRIPOS>PHARMACOPHORE\n", "@<TRIPOS>POINTS\n"])
        index = 0
        for element in pharmacophore:
            try:
                feat_name = pharmagist_element_name[element.feature_name]
            except KeyError:
                continue
            element_inx = str(index + 1)
            line = element_inx.rjust(7)
            line += " " + feat_name

            center = np.around(puw.get_value(element.center, to_unit="angstroms"), 4)
            x = _pad_coordinate_with_zeros(center[0]).rjust(16)
            y = _pad_coordinate_with_zeros(center[1]).rjust(10)
            z = _pad_coordinate_with_zeros(center[2]).rjust(10)

            line += x + y + z + " "
            line += pharmagist_element_specs[element.feature_name].rjust(5)
            line += str(index).rjust(5)
            line += pharmagist_element_specs[element.feature_name].rjust(6)
            line += "0.0000\n".rjust(12)
            lines.append(line)

            index += 1

        if save_props:
            props = _get_props_dict(pharmacophore)
            if len(props) > 0:
                lines.append("@<TRIPOS>PROPERTIES\n")
                for prop_name, value in props.items():
                    line = prop_name.rjust(len(prop_name) + 7) + " "
                    value = str(value)
                    line += (value + '\n').rjust(20 - len(prop_name))
                    lines.append(line)

    return lines


def save_mol2(pharmacophores, file_name, save_props=True):
    """ Save several pharmacophores to a mol2 file.

        Parameters
        ----------
        pharmacophores : list[Pharmacophore]
            A list of pharmacophores

        file_name : str
            Name of the file that will be created.

        save_props : bool
            Whether to save the properties stored in the pharmacophore.

    """
    contents = _mol2_pharmacophores(pharmacophores, save_props)
    with open(file_name, "w") as fp:
        fp.writelines(contents)

# PH4 (MOE) Files


def _ph4_pharmacophore(pharmacophoric_points):
    """ Returns a string with the necessary data to create a MOE ph4 file.

        Parameters
        ----------
        pharmacophoric_points : Pharmacophore
            List with the pharmacophoric points.


        Returns
        ------
        ph4_str : str
            The pharmacophore string.
    """
    oph_to_moe = {
        "aromatic ring": "Aro",
        "hydrophobicity": "Hyd",
        "hb acceptor": "Acc",
        "hb donor": "Don",
        "positive charge": "Cat",
        "negative charge": "Ani",
    }

    now = datetime.datetime.now()
    ph4_str = "#moe:ph4que" + " " + str(now.year) + "." + str(now.month) + "\n"
    ph4_str += "#pharmacophore 5 tag t value *\n"
    ph4_str += "scheme t Unified matchsize i 0 title t s $\n"
    ph4_str += f"#feature {len(pharmacophoric_points)} expr tt color ix x r y r z r r r ebits ix gbits ix\n"

    excluded_spheres = []
    for element in pharmacophoric_points:
        if element.feature_name == "excluded volume":
            excluded_spheres.append(element)
            continue
        feat_name = oph_to_moe[element.feature_name]
        center = puw.get_value(element.center, to_unit="angstroms")
        radius = puw.get_value(element.radius, to_unit="angstroms")
        radius_str = f"{radius:.3f}" + " "
        ph4_str += feat_name + " df2f2 "
        for ii in range(3):
            coord = f"{(center[ii]):.3f}" + " "
            ph4_str += coord
        ph4_str += radius_str + "0 300 "

    if excluded_spheres:
        ph4_str += "\n#volumesphere 90 x r y r z r r r\n"
        for excluded in excluded_spheres:
            center = puw.get_value(excluded.center, to_unit="angstroms")
            radius = puw.get_value(excluded.radius, to_unit="angstroms")
            radius_str = f"{radius:.3f}" + " "
            for ii in range(3):
                coord = f"{(center[ii]):.3f}" + " "
                ph4_str += coord
            ph4_str += radius_str

    ph4_str += "\n#endpharmacophore"

    return ph4_str


def save_ph4(pharmacophore, file_name):
    """ Save a pharmacophore to a ph4 (MOE) file.

        Parameters
        ----------
        pharmacophore : Pharmacophore
            A pharmacophore.

        file_name : str
            Name of the file that will be created.

    """
    contents = _ph4_pharmacophore(pharmacophore)
    with open(file_name, "w") as fp:
        fp.write(contents)

# JSON files


def _json_pharmacophore(pharmacophore):
    """ Returns a dictionary with the necessary data to construct a json pharmacophore file.

        Parameters
        ----------
        pharmacophore : Pharmacophore
            Pharmacophoric points that will be used to construct the dictionary

        Returns
        -------
        pharmacophore_dict : dict[str, Any]
            Dictionary with the necessary info to construct a .json pharmer file.
    """
    pharmer_element_name = {
        "aromatic ring": "Aromatic",
        "hydrophobicity": "Hydrophobic",
        "hb acceptor": "HydrogenAcceptor",
        "hb donor": "HydrogenDonor",
        "included volume": "InclusionSphere",
        "excluded volume": "ExclusionSphere",
        "positive charge": "PositiveIon",
        "negative charge": "NegativeIon",
    }
    points = []
    for element in pharmacophore:
        point_dict = {}
        temp_center = puw.get_value(element.center, to_unit='angstroms')
        point_dict["name"] = pharmer_element_name[element.feature_name]
        point_dict["svector"] = {}
        if element.has_direction:
            point_dict["hasvec"] = True
            point_dict["svector"]["x"] = element.direction[0]
            point_dict["svector"]["y"] = element.direction[1]
            point_dict["svector"]["z"] = element.direction[2]
        else:
            point_dict["hasvec"] = False
            point_dict["svector"]["x"] = 1.
            point_dict["svector"]["y"] = 0.
            point_dict["svector"]["z"] = 0.
        point_dict["x"] = temp_center[0]
        point_dict["y"] = temp_center[1]
        point_dict["z"] = temp_center[2]
        point_dict["radius"] = puw.get_value(element.radius, to_unit='angstroms')
        point_dict["enabled"] = True
        point_dict["vector_on"] = 0
        point_dict["minsize"] = ""
        point_dict["maxsize"] = ""
        point_dict["selected"] = False

        points.append(point_dict)

    return {"points": points}


def save_json(pharmacophore, file_name):
    """ Save a pharmacophore to a JSON file.

        Parameters
        ----------
        pharmacophore : Pharmacophore
            A pharmacophore.

        file_name : str
            Name of the file that will be created.

        """
    contents = _json_pharmacophore(pharmacophore)
    with open(file_name, "w") as fp:
        json.dump(contents, fp)

# PML (LigandScout) Files


def _set_coords_and_radius(tree_element, coords, radius, feat_id):
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


def _ligandscout_xml_tree(pharmacophore):
    """ Returns a xml element tree necessary to create a ligandscout pharmacophore file.

        Parameters
        ----------
        pharmacophore : openpharmacophore.Pharmacophore
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

    for ii, element in enumerate(pharmacophore):
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
            _set_coords_and_radius(point, (x, y, z),
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
            _set_coords_and_radius(plane, (x, y, z),
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
            _set_coords_and_radius(volume, (x, y, z),
                                  radius, feat_id)

    tree._setroot(document)

    return tree


def save_pml(pharmacophore, file_name):
    """ Save a pharmacophore to a pml (Ligand Scout) file.

        Parameters
        ----------
        pharmacophore : Pharmacophore
            A pharmacophore.

        file_name : str
            Name of the file that will be created.

        """
    xml_tree = _ligandscout_xml_tree(pharmacophore)
    xml_string = ET.tostring(xml_tree)
    with open(file_name, "wb") as fp:
        fp.write(xml_string)
