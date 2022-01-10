import pyunitwizard as puw
import numpy as np

def to_pharmagist(pharmacophores, file_name):
    """ Save pharmacophores to pharmagist mol2 format

        Parameters
        ----------
        pharmacophores : list of openpharmacophore.Pharmacophore or openpharmacophore.Pharmacophore
            Pharmacophore or pharmacophores that will be saved.

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
        pharmacophores : list of openpharmacophore.Pharmacophore or openpharmacophore.Pharmacophore
            Pharmacophore or pharmacophores that will be saved

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
        for i, element in enumerate(pharmacophore.elements):
            element_inx = str(i + 1)
            line += element_inx.rjust(7)
            try:
                feat_name = pharmagist_element_name[element.feature_name]
            except:
                continue
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
            line += str(i).rjust(5)
            line += pharmagist_element_specs[element.feature_name].rjust(6)
            line += "0.0000\n".rjust(12)
            lines.append(line)
            line = ""

        lines.append("@<TRIPOS>BOND\n")
        for l in lines:
            doc.append(l)
        
    return doc