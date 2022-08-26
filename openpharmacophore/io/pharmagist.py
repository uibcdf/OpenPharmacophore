from openpharmacophore import PharmacophoricPoint
from openpharmacophore.utils import bisection
import pyunitwizard as puw
import numpy as np
import re
from typing import List


def to_pharmagist(pharmacophores, file_name):
    """ Save pharmacophores to pharmagist mol2 format

        Parameters
        ----------
        pharmacophores : list of openpharmacophore.Pharmacophore or openpharmacophore.Pharmacophore
            The pharmacophores that will be saved.

        file_name : str
            Name of the file containing the pharmacophore.

        Notes
        -------
        Nothing is returned. A new file is written.
    """
    doc = _pharmagist_file_info(pharmacophores)
    with open(file_name, "w") as f:
        f.writelines(doc)


def _pharmagist_file_info(pharmacophores):
    """ Get necessary info to create a pharmagist mol2 file to store pharmacophores.

        Parameters
        ----------
        pharmacophores : list of Pharmacophore or Pharmacophore
            The pharmacophores that will be saved

        Returns
        -------
        doc : list of list
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

    doc = [] # list to store all pharmacophores
    for pharmacophore in pharmacophores:
        lines = ["@<TRIPOS>MOLECULE\n", "@<TRIPOS>ATOM\n"] # list to store all lines for a single pharmacophore
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
            if center[0] < 0:
                x = str(center[0]).ljust(7,"0").rjust(16)
            else:
                x = str(center[0]).ljust(6,"0").rjust(16)
            if center[1] < 0:
                y = str(center[1]).ljust(7,"0").rjust(10)
            else:
                y = str(center[1]).ljust(6,"0").rjust(10)
            if center[2] < 0:
                z = str(center[2]).ljust(7,"0").rjust(10)
            else:
                z = str(center[2]).ljust(6,"0").rjust(10)
            line += x + y + z + " "
            line += pharmagist_element_specs[element.feature_name].rjust(5)
            line += str(index).rjust(5)
            line += pharmagist_element_specs[element.feature_name].rjust(6)
            line += "0.0000\n".rjust(12)
            lines.append(line)
            line = ""
            
            index += 1

        lines.append("@<TRIPOS>BOND\n")
        for l in lines:
            doc.append(l)
        
    return doc


def load_pharmacophores(file_name: str) -> List[List[PharmacophoricPoint]]:
    """ Loads pharmacophores from a pharmagist mol2 file.

        Parameters
        ----------
        file_name : str
            Name of the file containing the pharmacophore.

        Returns
        -------
        pharmacophores : list of openpharmacophore.PharmacophoricPoint
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

    # This list will store a list with the pharmacophoric points for each pharmacophore
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
                bisection.insort_right(points, element, key=lambda p: p.short_name)
            if "@<TRIPOS>BOND" in line:
                pharmacophores.append(points)
    pharmacophores = [ph for ph in pharmacophores if len(ph) > 0]  # filter empty lists

    return pharmacophores
