from openpharmacophore import _puw, Pharmacophore
from openpharmacophore import pharmacophoric_elements as elements
import molsysmt as msm
import json
import pyunitwizard as puw


def from_pharmer(pharmacophore):

    tmp_pharmacophore = Pharmacophore()

    if type(pharmacophore)==str:
        if pharmacophore.endswith('.json'):
            with open(pharmacophore, "r") as fff:
                pharmacophore = json.load(fff)
        else:
            raise NotImplementedError

    def get_pharmer_element_properties(element, direction=False):
        center = _puw.quantity([element['x'], element['y'], element['z']], 'angstroms')
        radius = _puw.quantity(element['radius'], 'angstroms')
        if direction:
            direction = [element['svector']['x'], element['svector']['y'], element['svector']['z']]
            return center, radius, direction

        return center, radius

    for pharmer_element in pharmacophore['points']:
        pharmer_feature_name = pharmer_element['name']

        if pharmer_feature_name=='Aromatic':
            center, radius, direction = get_pharmer_element_properties(pharmer_element, direction=True)
            element = elements.AromaticRingSphereAndVector(center, radius, direction)
            tmp_pharmacophore.add_element(element)

        elif pharmer_feature_name=='Hydrophobic':
            center, radius = get_pharmer_element_properties(pharmer_element, direction=False)
            element = elements.HydrophobicSphere(center, radius)
            tmp_pharmacophore.add_element(element)

        elif pharmer_feature_name=='HydrogenAcceptor':
            center, radius, direction = get_pharmer_element_properties(pharmer_element, direction=True)
            element = elements.HBAcceptorSphereAndVector(center, radius, direction)
            tmp_pharmacophore.add_element(element)

        elif pharmer_feature_name=="HydrogenDonor":
            center, radius, direction = get_pharmer_element_properties(pharmer_element, direction=True)
            element = elements.HBDonorSphereAndVector(center, radius, direction)
            tmp_pharmacophore.add_element(element)

        elif pharmer_feature_name=="PositiveIon":
            center, radius = get_pharmer_element_properties(pharmer_element, direction=False)
            element = elements.PositiveChargeSphere(center, radius)
            tmp_pharmacophore.add_element(element)
        
        elif pharmer_feature_name=="NegativeIon":
            center, radius = get_pharmer_element_properties(pharmer_element, direction=False)
            element = elements.NegativeChargeSphere(center, radius)
            tmp_pharmacophore.add_element(element)

        elif pharmer_feature_name=="ExclusionSphere":
            center, radius = get_pharmer_element_properties(pharmer_element, direction=False)
            element = elements.ExcludedVolumeSphere(center, radius)
            tmp_pharmacophore.add_element(element)

        elif pharmer_feature_name=='InclusionSphere':
            center, radius = get_pharmer_element_properties(pharmer_element, direction=False)
            element = elements.IncludedVolumeSphere(center, radius)
            tmp_pharmacophore.add_element(element)

        else:
            raise NotImplementedError

    if "ligand" in pharmacophore:
        ligand = msm.convert(pharmacophore["ligand"], to_form="molsysmt.MolSys")
        #receptor = msm.convert(pharmacophore["receptor"], to_form="molsysmt.MolSys")
        #tmp_pharmacophore.molecular_system = msm.merge([ligand, receptor])

        tmp_pharmacophore.molecular_system = ligand

    return tmp_pharmacophore

def to_pharmer(pharmacophore, file_name):

    pharmer_element_name = { # dictionary to map openpharmacophore feature names to pharmer feature names
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
    for element in pharmacophore.elements:
        point_dict = {}
        temp_center = puw.get_value(element.center, to_unit='angstroms')
        point_dict["name"] = pharmer_element_name[element.feature_name]
        point_dict["svector"] = {}
        if hasattr(element, "direction"): 
            point_dict["hasvec"] = True
            point_dict["svector"]["x"] = element.direction[0]
            point_dict["svector"]["y"] = element.direction[1] 
            point_dict["svector"]["z"] = element.direction[2]  
        else: 
            point_dict["hasvec"] = False
            point_dict["svector"]["x"] = 1
            point_dict["svector"]["y"] = 0
            point_dict["svector"]["z"] = 0 
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

    pharmacophore = {}
    pharmacophore["points"] = points

    # TODO: add ligand and/or receptor
    
    with open(file_name, "w") as outfile:
        json.dump(pharmacophore, outfile)
