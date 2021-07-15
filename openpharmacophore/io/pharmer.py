from openpharmacophore import _puw, Pharmacophore
from openpharmacophore import pharmacophoric_elements as elements
import molsysmt as msm
import json

def from_pharmer(pharmacophore):

    tmp_pharmacophore = Pharmacophore()

    if type(pharmacophore)==str:
        if pharmacophore.endswith('.json'):
            with open(pharmacophore, "r") as fff:
                pharmacophore = json.load(fff)
        else:
            raise NotImplementedError

    for pharmer_element in pharmacophore['points']:
        pharmer_feature_name = pharmer_element['name']
        if pharmer_feature_name=='Aromatic':
            center = _puw.quantity([pharmer_element['x'], pharmer_element['y'], pharmer_element['z']], 'angstroms')
            radius = _puw.quantity(pharmer_element['radius'], 'angstroms')
            direction = [pharmer_element['svector']['x'], pharmer_element['svector']['y'], pharmer_element['svector']['z']]
            element = elements.AromaticRingSphereAndVector(center, radius, direction)
            tmp_pharmacophore.add_element(element)
        elif pharmer_feature_name=='Hydrophobic':
            center = _puw.quantity([pharmer_element['x'], pharmer_element['y'], pharmer_element['z']], 'angstroms')
            radius = _puw.quantity(pharmer_element['radius'], 'angstroms')
            element = elements.HydrophobicSphere(center, radius)
            tmp_pharmacophore.add_element(element)
        elif pharmer_feature_name=='HydrogenAcceptor':
            center = _puw.quantity([pharmer_element['x'], pharmer_element['y'], pharmer_element['z']], 'angstroms')
            radius = _puw.quantity(pharmer_element['radius'], 'angstroms')
            direction = [pharmer_element['svector']['x'], pharmer_element['svector']['y'], pharmer_element['svector']['z']]
            element = elements.HBAcceptorSphereAndVector(center, radius, direction)
            tmp_pharmacophore.add_element(element)
        elif pharmer_feature_name=='InclusionSphere':
            center = _puw.quantity([pharmer_element['x'], pharmer_element['y'], pharmer_element['z']], 'angstroms')
            radius = _puw.quantity(pharmer_element['radius'], 'angstroms')
            element = elements.IncludedVolumeSphere(center, radius)
            tmp_pharmacophore.add_element(element)
        else:
            raise NotImplementedError

    ligand = msm.convert(pharmacophore["ligand"], to_form="molsysmt.MolSys")
    #receptor = msm.convert(pharmacophore["receptor"], to_form="molsysmt.MolSys")
    #tmp_pharmacophore.molecular_system = msm.merge([ligand, receptor])

    tmp_pharmacophore.molecular_system = ligand

    return tmp_pharmacophore

def to_pharmer(pharmacophore, file_name=None):

    raise NotImplementedError

