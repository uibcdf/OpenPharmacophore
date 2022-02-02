from openpharmacophore import PharmacophoricPoint, LigandBasedPharmacophore
import pyunitwizard as puw
import re

def read_pharmagist(file_name, pharmacophore_index=None):
    """ Loads pharmacophores from a pharmagist mol2 file.

        Parameters
        ----------
        file_name : str
            Name of the file containing the pharmacophore.

        pharmacophore_index : int, optional
            If given, returns the pharmacophore which index is provided. 
            Else, all the pharmcophores from the file will be returned.

        Returns
        -------
        pharmacophores : list of openpharmacophore.LigandBasedPharmacophores
            A list of ligand based pharmacophores.
        
        pharmacophoric_points : list of openpharmacophore.PharmacophoricPoint
            A list of pharmacophoric points from a single pharmacophore.
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
                center = [float(coord) for coord in point_line[2: 5]] # convert coordinates to float
                center = puw.quantity(center, "angstroms")
                element = PharmacophoricPoint(
                    feat_type=feat_type,
                    center=center, 
                    radius=puw.quantity(1.0, "angstroms"))
                points.append(element)
            if "@<TRIPOS>BOND" in line:
                pharmacophores.append(points)
    pharmacophores = [ph for ph in pharmacophores if len(ph) > 0] # filter empty lists
    
    if pharmacophore_index is None:
        pharmacophores = [LigandBasedPharmacophore(ph) for ph in pharmacophores]
        return pharmacophores

    pharmacophoric_points = pharmacophores[pharmacophore_index]
    return pharmacophoric_points