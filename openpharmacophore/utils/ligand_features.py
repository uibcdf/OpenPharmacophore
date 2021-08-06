from rdkit import RDConfig, Chem
from rdkit.Chem import ChemicalFeatures
from openpharmacophore._private_tools.exceptions import NoConformersError, PointTypeError
from openpharmacophore.utils.centroid import feature_centroid
from openpharmacophore.utils.rdkit_to_point import rdkit_to_point
import numpy as np
import os

def rdkit_points(ligands, radius, feat_list=None, point_type="spheres"):
    """
        Get pharmacophoric points for a list of ligands using rdkit chemical feature definition. 

        Parameters
        ----------
        ligands: :obj: list of rdkit.Chem.rdchem.Mol
            List of ligands. Each ligand needs to have at least one conformer.
        
        radius: float
            Lenght of the radius of the parmacohporic points. Required if point type is 'spheres' or
            'spheres_vectors'.

        feat_list: list of str
            List of features that will be used to compute the pharmacophore.
        
        point_type: str
            Type of pharmacophoric points to be returned.

        Returns
        -------
        points: a list of openpharmacophore.pharmacophoric_elements

    """
    fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
    factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
    points = {}
        
    for i, ligand in enumerate(ligands):

        n_conformers = ligand.GetNumConformers()
        if n_conformers == 0:
            raise NoConformersError(n_conformers)

        ligand_id = "ligand_" + str(i)
        points[ligand_id] = {}

        feats = factory.GetFeaturesForMol(ligand)
        
        for f in feats:
            
            feat_name = f.GetFamily()
            if feat_name not in feat_list:
                continue
            atom_idxs = f.GetAtomIds()
            
            for conformer_idx in range(n_conformers):
                                    
                if len(atom_idxs) > 1: # Get the centroid of that feature
                    coords = feature_centroid(ligand, atom_idxs, conformer_idx)
                else:
                    position = ligand.GetConformer(conformer_idx).GetAtomPosition(atom_idxs[0])
                    coords = np.zeros((3,))
                    coords[0] = position.x
                    coords[1] = position.y
                    coords[2] = position.z

                # TODO: find direction vector for vectorial pharmacophoric points

                point = rdkit_to_point(feat_name, coords, radius=radius, point_type=point_type)
                conformer_id = "conformer_" + str(conformer_idx)

                if conformer_id not in points[ligand_id]:
                    points[ligand_id][conformer_id] = []
                points[ligand_id][conformer_id].append(point)
                
    return points

def custom_definition_points(ligands, radius, feat_list, feat_def, point_type="spheres"):
    """
        Get pharmacophoric points for a list of ligands using rdkit chemical feature definition. 

        Parameters
        ----------
        ligands: :obj: list of rdkit.Chem.rdchem.Mol
            List of ligands. Each ligand needs to have at least one conformer.
        
        radius: float
            Lenght of the radius of the parmacohporic points. Required if point type is 'spheres' or
            'spheres_vectors'.

        feat_list: list of str
            List of features that will be used to compute the pharmacophore.
        
        feat_def: dict
            Definitions of the pharmacophoric points. 
            Dictionary which keys are SMARTS strings and values are feature names.
        
        point_type: str
            Type of pharmacophoric points to be returned.

        Returns
        -------
        points: a list of openpharmacophore.pharmacophoric_elements

    """
    
    points = {}

    for i, ligand in enumerate(ligands):
        n_conformers = ligand.GetNumConformers()
        if n_conformers == 0:
            raise NoConformersError(n_conformers)
        
        ligand_id = "ligand_" + str(i)
        points[ligand_id] = {}

        for feat, feat_type in feat_def.items():
            if feat_type not in feat_list:
                continue
            pattern = Chem.MolFromSmarts(feat)  
            atom_idxs = ligand.GetSubstructMatch(pattern)
            if len(atom_idxs) == 0:
                continue
            for conformer_idx in range(n_conformers):
                if len(atom_idxs) == 1: # Donor or acceptor feature
                    position = ligand.GetConformer(conformer_idx).GetAtomPosition(atom_idxs[0])
                    coords = np.zeros((3,))
                    coords[0] = position.x
                    coords[1] = position.y
                    coords[2] = position.z
                else: # Aromatic, hydrophobic, positive or negative feature
                    coords = feature_centroid(ligand, atom_idxs, conformer_idx)

                point = rdkit_to_point(feat_type, coords, radius=radius, point_type=point_type)
                conformer_id = "conformer_" + str(conformer_idx)

                if conformer_id not in points[ligand_id]:
                    points[ligand_id][conformer_id] = []
                points[ligand_id][conformer_id].append(point)

    return points


def ligands_pharmacophoric_points(ligands, radius, feat_list=None, feat_def=None, 
                                point_type="spheres"):

    """
        Get pharmacophoric points for each ligand in a list of ligands. If a ligand has 
        more than one conformer, pharamcophoric points will be computed for each one.  

        Parameters
        ----------
        ligands: :obj: rdkit.Chem.rdmolfiles.SmilesMolSupplier or list of rdkit.Chem.rdchem.Mol
            List of ligands. Each ligand needs to have at least one conformer.
        
        radius: float
            Lenght of the radius of the parmacohporic points. Required if point type is 'spheres' or
            'spheres_vectors'.

        feat_list: list of str (optional)
            List of features that will be used to compute the pharmacophore.

        feat_def: dict (optional)
            Definitions of the pharmacophoric points. 
            Dictionary which keys are SMARTS strings and values are feature names.
        
        point_type: str
            Type of pharmacophoric points to be returned.

        Returns
        -------
        points: a list of openpharmacophore.pharmacophoric_elements

    """

    point_type_list = ["spheres", "spheres_vectors", "gaussian", "shapelet"]
    if point_type not in point_type_list:
        raise PointTypeError(f"Invalid point type. \"{point_type}\" is not a valid point type")

    if (point_type == "spheres" or point_type == "spheres_vectors") and radius is None:
        raise ValueError("Radius cannot be null if point type is spheres or spheres_vectors") 

    if isinstance(ligands, Chem.rdchem.Mol): # Check if it's a single molecule
        ligands = [ligands]
    elif not isinstance(ligands, list):
        ligands = list(ligands)
    
    if not feat_list:
        feat_list = ['Acceptor', 'Aromatic', 'Donor', 'Hydrophobe', 'PosIonizable', 'NegIonizable']
    
    if feat_def is None: # If no feature definition is given use rdkit one
        points = rdkit_points(ligands=ligands, radius=radius, feat_list=feat_list, point_type=point_type)
    else:
        points = custom_definition_points(ligands=ligands, radius=radius, feat_list=feat_list, feat_def=feat_def, point_type=point_type)
    
    return points

    
