import json
import re
import xml.etree.ElementTree as ET
import pyunitwizard as puw
from rdkit import Chem

from openpharmacophore import PharmacophoricPoint


# MOL2 Files


def read_mol2(file_name):
    """ Loads pharmacophores from a pharmagist mol2 file.

        Parameters
        ----------
        file_name : str
            Name of the file containing the pharmacophore.

        Returns
        -------
        pharmacophores : list[list[PharmacophoricPoint]]
            A list of ligand based pharmacophores.

    """
    # dictionary to map pharmagist feature names to openpharmacophore pharmacophoric_points
    pharmagist_element_name = {
        "AR": "aromatic ring",
        "HYD": "hydrophobicity",
        "ACC": "hb acceptor",
        "DON": "hb donor",
        "CAT": "positive charge",
        "ANI": "negative charge",
    }

    # This list will store a list with the pharmacophoric points of each pharmacophore
    pharmacophores = []
    pattern = r"[A-Z]{2,3} "
    with open(file_name, "r") as f:
        for line in f.readlines():
            match = re.search(pattern, line)
            if "@<TRIPOS>ATOM" in line:
                points = []
            if match:
                point_line = [p for p in line.split(" ") if p != ""]
                feat_type = pharmagist_element_name[point_line[1]]
                center = [float(coord) for coord in point_line[2: 5]]  # convert coordinates to float
                center = puw.quantity(center, "angstroms")
                element = PharmacophoricPoint(
                    feat_type=feat_type,
                    center=center,
                    radius=puw.quantity(1.0, "angstroms"))
                points.append(element)
            if "@<TRIPOS>BOND" in line:
                pharmacophores.append(points)
    pharmacophores = [ph for ph in pharmacophores if len(ph) > 0]  # filter empty lists

    return pharmacophores


# PH4 Files

def read_ph4(file_name: str):
    """ Loads a list of pharmacophoric points from a MOE ph4 file.

        Parameters
        ----------
        file_name : str
            Name of the file containing the pharmacophore.

        Returns
        -------
        points : list[PharmacophoricPoint]
            A list of pharmacophoric points.
    """
    moe_to_oph = {
        "Cat": "positive charge",
        "Ani": "negative charge",
        "Don": "hb donor",
        "Acc": "hb acceptor",
        "Hyd": "hydrophobicity",
        "Aro": "aromatic ring",
    }

    points = []

    with open(file_name, "r") as f:
        moe_ph4_string = f.read()

    # pieces = re.split("\s", moe_ph4_string)
    pieces = moe_ph4_string.split()
    pattern = r"Don|Acc|Aro|Cat|Ani|Hyd"
    # Index where the excluded volumes spheres coordinates, if any, start
    exclusion_inx = None
    for i, piece in enumerate(pieces):
        if re.match(pattern, piece):
            if len(piece) == 7:
                # Double features have seven characters i.e. Acc|Aro
                feat_name_1 = piece[0: 3]
                feat_name_2 = piece[4:]
                # Coordinates are always two items after the feature name
                x = float(pieces[i + 2])
                y = float(pieces[i + 3])
                z = float(pieces[i + 4])
                radius = float(pieces[i + 5])
                point_1 = PharmacophoricPoint(
                    feat_type=moe_to_oph[feat_name_1],
                    center=puw.quantity([x, y, z], "angstroms"),
                    radius=puw.quantity(radius, "angstroms")
                )
                point_2 = PharmacophoricPoint(
                    feat_type=moe_to_oph[feat_name_2],
                    center=puw.quantity([x, y, z], "angstroms"),
                    radius=puw.quantity(radius, "angstroms")
                )
                points.append(point_1)
                points.append(point_2)
            else:
                # Regular features and
                # features with weird names i.e. Cat$mDon, Acc2
                feat_name = piece[0:3]
                # Coordinates are always two items after the feature name
                x = float(pieces[i + 2])
                y = float(pieces[i + 3])
                z = float(pieces[i + 4])
                radius = float(pieces[i + 5])
                point = PharmacophoricPoint(
                    feat_type=moe_to_oph[feat_name],
                    center=puw.quantity([x, y, z], "angstroms"),
                    radius=puw.quantity(radius, "angstroms")
                )
                points.append(point)

        if piece == "#volumesphere":
            # Excluded volume spheres are 10 items after #volumesphere
            exclusion_inx = i + 10
            break

    if exclusion_inx:
        count = 1
        for p in pieces[exclusion_inx:]:
            if p == "#volume":
                break
            if count == 1:
                x = float(p)
                count += 1
            elif count == 2:
                y = float(p)
                count += 1
            elif count == 3:
                z = float(p)
                count += 1
            elif count == 4:
                radius = float(p)
                excluded_sphere = PharmacophoricPoint(
                    feat_type="excluded volume",
                    center=puw.quantity([x, y, z], "angstroms"),
                    radius=puw.quantity(radius, "angstroms")
                )
                points.append(excluded_sphere)
                count = 1

    return points


# JSON Files

def _get_pharmer_element_properties(element, direction=False):
    center = puw.quantity([element['x'], element['y'], element['z']], 'angstroms')
    radius = puw.quantity(element['radius'], 'angstroms')
    if direction:
        direction = [element['svector']['x'], element['svector']['y'], element['svector']['z']]
        return center, radius, direction

    return center, radius


def read_json(pharmacophore_file, load_mol_sys=False):
    """ Loads a pharmacophore from a json file

        Parameters
        ----------
        pharmacophore_file : str
            name of the file containing the pharmacophore

        load_mol_sys : bool
            If true loads the molecular system associated to the pharmacophore (Default: False).

        Returns
        -------
        points : list[PharmacophoricPoint]
            A list of pharmacophoric points.

        molecular_system : rdkit.Mol
            The molecular system associated with the pharmacophore. If there is no molecular system or
            if load_mol_sys is set to false, None is returned.

        ligand : rdkit.Mol
            The ligand associated to the pharmacophore in case there is one.
             If there is no ligand None is returned.
    """
    points = []
    molecular_system = None
    ligand = None

    if pharmacophore_file.endswith('.json'):
        with open(pharmacophore_file, "r") as fff:
            pharmacophore = json.load(fff)
    else:
        raise NotImplementedError

    feature_name = {
        "Aromatic": "aromatic ring",
        "Hydrophobic": "hydrophobicity",
        "HydrogenAcceptor": "hb acceptor",
        "HydrogenDonor": "hb donor",
        "InclusionSphere": "included volume",
        "ExclusionSphere": "excluded volume",
        "PositiveIon": "positive charge",
        "NegativeIon": "negative charge",
    }

    direction_points = ["Aromatic", "Hydrophobic", "HydrogenAcceptor", "HydrogenDonor"]

    for pharmer_element in pharmacophore['points']:

        pharmer_feature_name = pharmer_element['name']
        if pharmer_feature_name in direction_points:
            center, radius, direction = _get_pharmer_element_properties(pharmer_element, direction=True)
            element = PharmacophoricPoint(feature_name[pharmer_feature_name], center, radius, direction)
        else:
            center, radius = _get_pharmer_element_properties(pharmer_element, direction=False)
            element = PharmacophoricPoint(feature_name[pharmer_feature_name], center, radius)

        points.append(element)

    if load_mol_sys:
        try:
            if pharmacophore["ligand"]:
                ligand = Chem.MolFromPDBBlock(pharmacophore["ligand"])
        except KeyError:
            pass

        try:
            if pharmacophore["receptor"]:
                molecular_system = Chem.MolFromPDBBlock(pharmacophore["receptor"])
        except KeyError:
            pass

    return points, molecular_system, ligand


# PML (LigandScout) Files


_ligandscout_to_oph = {
    "PI": "positive charge",
    "NI": "negative charge",
    "HBD": "hb donor",
    "HBA": "hb acceptor",
    "H": "hydrophobicity",
    "AR": "aromatic ring",
    "exclusion": "excluded volume"
}


def _pharmacophoric_point_from_point_or_volume(element):
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

    return PharmacophoricPoint(_ligandscout_to_oph[feat_name],
                               puw.quantity([x, y, z], "angstroms"),
                               puw.quantity(radius, "angstroms"))


def _pharmacophoric_point_from_vector_or_plane(element):
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
    return PharmacophoricPoint(_ligandscout_to_oph[feat_name],
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
            point = _pharmacophoric_point_from_point_or_volume(element)

        elif element.tag == "vector" or element.tag == "plane":
            point = _pharmacophoric_point_from_vector_or_plane(element)

        else:
            continue

        points.append(point)

    return points
