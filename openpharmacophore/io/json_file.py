from ..pharmacophore import PharmacophoricPoint
from rdkit import Chem
import pyunitwizard as puw
import json


def get_pharmer_element_properties(element, direction=False):
    center = puw.quantity([element['x'], element['y'], element['z']], 'angstroms')
    radius = puw.quantity(element['radius'], 'angstroms')
    if direction:
        direction = [element['svector']['x'], element['svector']['y'], element['svector']['z']]
        return center, radius, direction

    return center, radius


def load_json_pharmacophore(pharmacophore_file, load_mol_sys=False):
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
            center, radius, direction = get_pharmer_element_properties(pharmer_element, direction=True)
            element = PharmacophoricPoint(feature_name[pharmer_feature_name], center, radius, direction)
        else:
            center, radius = get_pharmer_element_properties(pharmer_element, direction=False)
            element = PharmacophoricPoint(feature_name[pharmer_feature_name], center, radius)

        points.append(element)

    if load_mol_sys:
        try:
            if pharmacophore["ligand"]:
                ligand = Chem.rdmolfiles.MolFromPDBBlock(pharmacophore["ligand"])
        except KeyError:
            pass

        try:
            if pharmacophore["receptor"]:
                molecular_system = Chem.rdmolfiles.MolFromPDBBlock(pharmacophore["receptor"])
        except KeyError:
            pass

    return points, molecular_system, ligand


def json_pharmacophoric_elements(pharmacophoric_points):
    """ Returns a dictionary with the necessary data to construct a json pharmacophore file.

        Parameters
        ----------
        pharmacophoric_points : list[PharmacophoricPoint]
            Pharmacophore points that will be used to construct the dictionary

        Returns
        -------
        pharmacophore_dict : dict
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
    for element in pharmacophoric_points:
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
