from ..pharmacophore import PharmacophoricPoint
import pyunitwizard as puw
import numpy as np
import re


def pad_coordinate_with_zeros(coord):
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


def mol2_file_info(pharmacophores):
    """ Get necessary info to create a mol2 file to store pharmacophores.

        Parameters
        ----------
        pharmacophores : list[list[PharmacophoricPoint]]
            Nested list of pharmacophoric points were each entry represents a pharmacophore

        Returns
        -------
        doc : list[list[str]]
            List where each sublist contains a pharmacophore represented as a mol2 string.
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

    if not isinstance(pharmacophores, list):
        pharmacophores = [pharmacophores]

    doc = []  # list to store all pharmacophores
    for pharmacophore in pharmacophores:
        lines = ["@<TRIPOS>MOLECULE\n", "@<TRIPOS>ATOM\n"]  # list to store all lines for a single pharmacophore
        line = ""
        index = 0
        for element in pharmacophore:
            try:
                feat_name = pharmagist_element_name[element.feature_name]
            except KeyError:
                continue
            element_inx = str(index + 1)
            line += element_inx.rjust(7)
            line += " " + feat_name
            # Get point coordinates
            center = np.around(puw.get_value(element.center, to_unit="angstroms"), 4)
            # Pad coordinates with zeros to the right. Number of zeros depends on sign
            x = pad_coordinate_with_zeros(center[0]).rjust(16)
            y = pad_coordinate_with_zeros(center[1]).rjust(10)
            z = pad_coordinate_with_zeros(center[2]).rjust(10)
            line += x + y + z + " "
            line += pharmagist_element_specs[element.feature_name].rjust(5)
            line += str(index).rjust(5)
            line += pharmagist_element_specs[element.feature_name].rjust(6)
            line += "0.0000\n".rjust(12)
            lines.append(line)
            line = ""

            index += 1

        lines.append("@<TRIPOS>BOND\n")
        for line in lines:
            doc.append(line)

    return doc


def load_mol2_pharmacophoric_points(file_name):
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
